#!/usr/bin/env python3
"""Stable request fingerprints. Files use canonical names and nanosecond stat data."""
import argparse
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--file', action='append', default=[])
p.add_argument('--text-file', action='append', default=[])
p.add_argument('--value', action='append', default=[])
a = p.parse_args()
files = []
for name in a.file + a.text_file:
    path = Path(name).resolve(strict=True)
    s = path.stat()
    row = [str(path), s.st_size, s.st_mtime_ns]
    if name in a.text_file:
        row.append(hashlib.sha256(path.read_bytes()).hexdigest())
    files.append(row)
print(json.dumps({'schema': 2, 'values': a.value, 'files': files}, sort_keys=True))
