#!/usr/bin/env python3
"""Compact terminal output and retain the complete workload log in the workflow logs directory.

All workflow log styles share this launcher.
It runs inside the same verified cgroup as the workload and all descendants.
"""
from __future__ import annotations

import argparse
from collections import deque
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
    with log_path.open("ab", buffering=0) as log:
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


def log_label(value):
    """Accept workflow names without allowing paths in the log filename."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise argparse.ArgumentTypeError("label must start with a letter/digit and contain only letters, digits, _, . or -")
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--label", type=log_label, required=True, help="Workflow label used in the log filename")
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
        path = a.log_dir / f"{a.label}.log"
        return run(a.command, path, a.interval, a.verbose)
    except OSError as exc:
        print(f"allele3 console: {exc}", file=sys.stderr)
        return 2



import json
class HervScreenLog(ScreenLog):
    def __init__(self, interval, verbose=False):
        super().__init__(interval, verbose)
        self.tail = deque(maxlen=3)
        self.prepared = set()
        self.annotated = 0
        self.total_chunks = None
        self.warnings = set()

    def progress(self, stage, message, complete=False):
        now = time.monotonic()
        if complete or now - self.last.get(stage, float('-inf')) >= self.interval:
            print(message, flush=True)
            self.last[stage] = now

    def line(self, raw):
        if self.verbose:
            return  # feed() already writes the unmodified bytes.
        line = raw.decode('utf-8', errors='replace').strip()
        if not line:
            return
        if line.startswith(('[dataset]', '[build-check]')):
            print(line, flush=True)
            return
        if line.startswith('{'):
            try:
                event = json.loads(line)
            except (ValueError, TypeError):
                return
            if not isinstance(event, dict):
                return
            if event.get('event') in ('prepared', 'prepared_pgen'):
                before = len(self.prepared)
                self.prepared.add(event.get('chromosome'))
                complete = before < 22 and len(self.prepared) == 22
                self.progress('prepare', f'[herv] prepared {len(self.prepared)}/22 chromosomes; latest={event.get("chromosome")}', complete)
            elif event.get('event') == 'annotated':
                self.annotated += 1
                if self.total_chunks is None and os.environ.get('HERV_SCREEN_OUT'):
                    try:
                        self.total_chunks = json.loads((Path(os.environ['HERV_SCREEN_OUT']) / 'plan.json').read_text())['chunks']
                    except (OSError, ValueError, KeyError):
                        pass
                count = f'{self.annotated}/{self.total_chunks}' if self.total_chunks else str(self.annotated)
                self.progress('annotate', f'[herv] annotated {count} chunks; latest={event.get("chromosome")} {event.get("label")}', self.annotated == self.total_chunks)
            return
        if line.startswith('herv:') or re.search(r'\b(error|fatal|failed|fail|traceback)\b', line, re.I):
            if 'Command ' in line and 'returned non-zero exit status' in line:
                status = re.search(r'exit status (\d+)', line)
                program = 'RepeatMasker' if "'/usr/bin/time'" in line else 'HERV stage'
                line = f'{program} failed (exit {status[1] if status else "unknown"}); details in the full log'
            if line not in self.tail:
                self.tail.append(line)
            return
        if re.search(r'\bwarning\b', line, re.I):
            if line not in self.warnings:
                self.warnings.add(line)
                if len(self.warnings) <= 3:
                    print(line, flush=True)
            return
        if re.match(r'^\[herv\] (START|DONE|VERIFY|REUSE|PREFLIGHT)\b', line):
            print(line, flush=True)



class PipelineScreenLog(ScreenLog):
    def __init__(self, interval, verbose=False):
        super().__init__(interval, verbose)
        self.tail = deque(maxlen=5)
        self.diagnostics = set()

    def line(self, raw):
        if self.verbose:
            return  # feed() already forwards unmodified output.
        line = raw.decode("utf-8", errors="replace").strip()
        if not line:
            return
        if line.startswith(("[dataset]", "[build-check]")):
            print(line, flush=True)
            return
        if line.startswith("[stage] "):
            stage = line.split(":", 1)[0]
            print(f"{stage}: running", flush=True)
            self.last.clear()
            self.diagnostics.clear()
            return
        if line.startswith(("[reuse] ", "[not_available] ", "[herv-focus] ")):
            print(line, flush=True)
            return
        # JSON, native banners and tracebacks are retained in the full log.
        # Python exception endings and our CLI error prefixes provide context.
        if line.startswith(("{", "}", '"')):
            return
        error = re.match(r"^(?:herv-[\w-]+:|[\w.]+(?:Error|Exception):)", line)
        if error or DIAGNOSTIC.search(line):
            summary = line[:500]
            if summary not in self.tail:
                self.tail.append(summary)
            if summary not in self.diagnostics and len(self.diagnostics) < 3:
                self.diagnostics.add(summary)
                print(summary, flush=True)
            return
        if line.startswith("["):
            now = time.monotonic()
            if now - self.last.get("progress", float("-inf")) >= self.interval:
                print(line[:500], flush=True)
                self.last["progress"] = now



if __name__ == '__main__':
    style = os.environ.get('ALLELE_LOG_STYLE', 'default')
    ScreenLog = {'herv': HervScreenLog, 'pipeline': PipelineScreenLog}.get(style, ScreenLog)
    raise SystemExit(main())
