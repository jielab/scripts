#!/usr/bin/env python3
"""HERV progress summaries; complete raw output remains in the workload log."""
from collections import deque
import importlib.util
import json
import os
from pathlib import Path
import re
import time

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('workload_console', ROOT / 'f/0.log.py')
console = importlib.util.module_from_spec(spec)
spec.loader.exec_module(console)


class HervScreenLog(console.ScreenLog):
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


if __name__ == '__main__':
    console.ScreenLog = HervScreenLog
    raise SystemExit(console.main())
