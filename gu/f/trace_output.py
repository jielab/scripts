#!/usr/bin/env python3
"""Publish TRACE NPZ outputs only after successful execution and ZIP validation."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import zipfile
import zlib


def valid_npz(path):
    try:
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            return bool(names) and all(n.endswith('.npy') for n in names) and archive.testzip() is None
    except (OSError, ValueError, EOFError, RuntimeError, zipfile.BadZipFile, zlib.error):
        return False


def signature(paths):
    return [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in paths]


def reusable(paths, receipt):
    try:
        if json.loads(receipt.read_text()) == signature(paths):
            return True
    except (OSError, ValueError):
        pass
    # Adopt complete pre-existing outputs after validating every ZIP member.
    if not all(valid_npz(p) for p in paths):
        return False
    seal(paths, receipt)
    return True


def seal(paths, receipt):
    part = receipt.with_name(receipt.name + '.part')
    part.write_text(json.dumps(signature(paths)) + '\n')
    part.replace(receipt)


def run(prefix, suffixes, command):
    prefix = Path(prefix)
    paths = [Path(str(prefix) + s) for s in suffixes]
    receipt = Path(str(prefix) + '.complete.json')
    if reusable(paths, receipt):
        print(f'[GU TRACE] SKIP verified={prefix}', flush=True)
        return 0
    receipt.unlink(missing_ok=True)
    # Delete invalid legacy files so a failed replacement cannot be adopted.
    for p in paths:
        p.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix='.trace-part-', dir=prefix.parent) as tmp:
        staged = Path(tmp) / prefix.name
        print(f'[GU TRACE] START output={prefix}', flush=True)
        proc = subprocess.Popen(command + ['-o', str(staged)], stdin=subprocess.DEVNULL,
                                start_new_session=True)
        try:
            rc = proc.wait()
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
        if rc:
            return rc if rc > 0 else 128 - rc
        outputs = [Path(str(staged) + s) for s in suffixes]
        if not all(valid_npz(p) for p in outputs):
            print(f'ERROR: incomplete TRACE NPZ output: {prefix}', flush=True)
            return 1
        for source, target in zip(outputs, paths):
            source.replace(target)
        seal(paths, receipt)
        print(f'[GU TRACE] DONE output={prefix}', flush=True)
    return 0


def interrupted(signum, frame):
    raise KeyboardInterrupt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--suffix', action='append', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('a TRACE command is required after --')
    signal.signal(signal.SIGTERM, interrupted)
    try:
        raise SystemExit(run(args.prefix, args.suffix, command))
    except KeyboardInterrupt:
        raise SystemExit(130)
