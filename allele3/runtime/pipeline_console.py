#!/usr/bin/env python3
"""Compact allele3 analysis output; raw bytes stay in the full workload log."""
from collections import deque
import importlib.util
from pathlib import Path
import re
import time

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("workload_console", ROOT / "f/0.log.py")
console = importlib.util.module_from_spec(spec)
spec.loader.exec_module(console)


class PipelineScreenLog(console.ScreenLog):
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
        if error or console.DIAGNOSTIC.search(line):
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


if __name__ == "__main__":
    console.ScreenLog = PipelineScreenLog
    raise SystemExit(console.main())
