#!/usr/bin/env python3
"""Compact terminal output and retain the complete workload log outside f/.

Keeping this launcher outside the analysis package preserves scan cache hashes.
It runs inside the same verified cgroup as the workload and all descendants.
"""
from __future__ import annotations

import argparse
from collections import deque
from datetime import datetime
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time


PROGRESS = re.compile(
    r"^\[(scan [^]]+|map|selection|selection index|modern intersect|PVAR|PGEN cache)\]"
)
DIAGNOSTIC = re.compile(r"\b(\w*warnings?|error|fatal|failed|traceback)\b", re.IGNORECASE)


class ScreenLog:
    def __init__(self, interval, verbose=False):
        self.interval = interval
        self.verbose = verbose
        self.last = {}
        self.pending = b""
        self.tail = deque(maxlen=30)

    def line(self, raw):
        line = raw.decode("utf-8", errors="replace").strip()
        if not line:
            return
        self.tail.append(line[-4096:])
        if self.verbose:
            return
        if DIAGNOSTIC.search(line):
            print(line, flush=True)
            return
        progress = PROGRESS.match(line)
        if progress:
            # Concurrent chromosomes share one screen budget, not one each.
            key = progress[1]
            if key.startswith("scan "):
                key = "scan"
            now = time.monotonic()
            if now - self.last.get(key, float("-inf")) < self.interval:
                return
            self.last[key] = now
        elif not line.startswith("["):
            # Native tool banners, PLINK progress and JSON remain in the file.
            return
        print(line, flush=True)

    def feed(self, chunk):
        if self.verbose:
            sys.stdout.buffer.write(chunk)
            sys.stdout.buffer.flush()
        parts = re.split(rb"[\r\n]", self.pending + chunk)
        self.pending = parts.pop()
        for part in parts:
            self.line(part)
        # Bound memory even for tools producing very long unterminated lines.
        if len(self.pending) >= 65536:
            self.line(self.pending)
            self.pending = b""

    def finish(self):
        if self.pending:
            self.line(self.pending)
            self.pending = b""


def run(command, log_path, interval=300, verbose=False):
    label = os.environ.get("WORKFLOW_NAME", "allele3")
    screen = ScreenLog(interval, verbose)
    started = time.monotonic()
    with log_path.open("xb", buffering=0) as log:
        print(f"[{label}] full log={log_path}", flush=True)
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            env=dict(os.environ, PYTHONUNBUFFERED="1"), start_new_session=True,
        )

        def forward(signum, _frame):
            try:
                os.killpg(process.pid, signum)
            except ProcessLookupError:
                pass

        handlers = {s: signal.signal(s, forward) for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}
        try:
            with process.stdout:
                while chunk := os.read(process.stdout.fileno(), 65536):
                    log.write(chunk)
                    screen.feed(chunk)
            code = process.wait()
            screen.finish()
        finally:
            if process.poll() is None:
                forward(signal.SIGTERM, None)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    forward(signal.SIGKILL, None)
                    process.wait()
            for signum, handler in handlers.items():
                signal.signal(signum, handler)
    if code:
        if not verbose:
            print(f"[{label}] last workload output:", file=sys.stderr)
            for line in screen.tail:
                print(line, file=sys.stderr)
        print(f"[{label}] failed exit={code}; full log={log_path}", file=sys.stderr, flush=True)
    else:
        print(f"[{label}] complete; elapsed={(time.monotonic() - started) / 60:.1f} min; full log={log_path}", flush=True)
    return code if code >= 0 else 128 - code


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--label", choices=("pipeline", "run", "scan", "evolve", "herv", "pca3d"), required=True)
    parser.add_argument("--interval", type=int, default=300)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    a = parser.parse_args(argv)
    if a.command[:1] == ["--"]:
        a.command.pop(0)
    if not a.command or a.interval <= 0:
        parser.error("a command and a positive screen interval are required")
    try:
        a.log_dir.mkdir(parents=True, exist_ok=True)
        path = a.log_dir / f"{a.label}_{datetime.now():%Y%m%d_%H%M%S_%f}_{os.getpid()}.log"
        return run(a.command, path, a.interval, a.verbose)
    except OSError as exc:
        print(f"allele3 console: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
