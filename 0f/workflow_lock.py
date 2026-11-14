"""Whole-file Linux OFD record locks, shared across GPFS compute nodes.

QiMing's BSD flock is node-local. OFD record locks were verified across two
nodes on its /work GPFS mount. Unlike process-scoped lockf, closing an unrelated
descriptor cannot release these locks, and threads with separate opens contend.
Never unlink a lock file: all readers/writers must retain the same inode.
Unsupported filesystems fail closed; there is no unsafe flock fallback.
"""
from contextlib import contextmanager
import argparse
import ctypes
import errno
import fcntl
import os
from pathlib import Path
import sys
import time


class _Flock(ctypes.Structure):
    _fields_ = [("type", ctypes.c_short), ("whence", ctypes.c_short),
                ("start", ctypes.c_longlong), ("length", ctypes.c_longlong),
                ("pid", ctypes.c_int)]


def flock(file, operation):
    """fcntl.flock-compatible flags using a descriptor-scoped record lock."""
    fd = file if isinstance(file, int) else file.fileno()
    mode = operation & ~fcntl.LOCK_NB
    if mode not in (fcntl.LOCK_SH, fcntl.LOCK_EX, fcntl.LOCK_UN):
        raise ValueError("Invalid file-lock operation")
    # Conda Python may have been compiled against pre-3.15 Linux headers even
    # when the running kernel supports OFD locks. 37 is the Linux LP64 UAPI
    # command (also exercised by the cross-node deployment probe).
    command = getattr(fcntl, "F_OFD_SETLK", 37 if sys.platform == 'linux'
                      and ctypes.sizeof(ctypes.c_long) == 8 else None)
    if command is None:
        raise RuntimeError("Linux OFD record locks are required for cross-node cache safety")
    kind = {fcntl.LOCK_SH: fcntl.F_RDLCK, fcntl.LOCK_EX: fcntl.F_WRLCK,
            fcntl.LOCK_UN: fcntl.F_UNLCK}[mode]
    record = bytes(_Flock(kind, os.SEEK_SET, 0, 0, 0))
    # Poll nonblocking: bounded waits allow Python signal handlers to run, and
    # long lock waits remain visible in scheduler logs.
    started = last = time.monotonic()
    while True:
        try:
            fcntl.fcntl(fd, command, record)
            return
        except OSError as error:
            if error.errno not in (errno.EAGAIN, errno.EACCES):
                raise
            if operation & fcntl.LOCK_NB:
                raise BlockingIOError(errno.EAGAIN, "Shared filesystem lock is held") from error
            now = time.monotonic()
            if now - last >= 60:
                try:
                    name = os.readlink(f"/proc/self/fd/{fd}")
                except OSError:
                    name = str(fd)
                print(f"[lock] waiting {now-started:.0f}s: {name}", flush=True)
                last = now
            time.sleep(.25)


@contextmanager
def file_lock(path, *, exclusive=True, wait=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Read locks also use O_RDWR: the same descriptor may later be upgraded.
    with path.open("a+") as stream:
        flock(stream, (fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
              | (0 if wait else fcntl.LOCK_NB))
        yield stream


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--shared", action="store_true")
    parser.add_argument("-x", "--exclusive", action="store_true")
    parser.add_argument("-u", "--unlock", action="store_true")
    parser.add_argument("-n", "--nonblocking", action="store_true")
    parser.add_argument("fd", type=int, help="Inherited open Bash descriptor (open with <>)")
    args = parser.parse_args()
    mode = fcntl.LOCK_UN if args.unlock else fcntl.LOCK_SH if args.shared else fcntl.LOCK_EX
    try:
        flock(args.fd, mode | (fcntl.LOCK_NB if args.nonblocking else 0))
    except BlockingIOError:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
