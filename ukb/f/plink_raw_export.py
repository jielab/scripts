#!/usr/bin/env python3
"""Stream PLINK's sample-major .raw output into gzip without a plain disk copy."""
import argparse
import gzip
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile


def export(output, command):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.raw-export.', dir=output.parent) as temporary:
        prefix = Path(temporary)/'data'
        fifo = Path(str(prefix)+'.raw')
        os.mkfifo(fifo)
        # Hold the pipe open until PLINK exits, including failure before its first write.
        keeper = os.open(fifo, os.O_RDWR)
        reader = fifo.open('rb')
        packed = Path(temporary)/'raw.gz'
        gzip_proc = plink_proc = None
        try:
            with packed.open('wb') as dest:
                gzip_proc = subprocess.Popen(['gzip', '-1c'], stdin=reader, stdout=dest)
                plink_proc = subprocess.Popen([*command, '--out', str(prefix)])
                while plink_proc.poll() is None:
                    if gzip_proc.poll() is not None:
                        raise RuntimeError('gzip exited before PLINK completed')
                    try: plink_proc.wait(timeout=0.2)
                    except subprocess.TimeoutExpired: pass
                rc = plink_proc.returncode
                os.close(keeper); keeper = None
                if rc: raise subprocess.CalledProcessError(rc, command)
                if gzip_proc.wait(): raise RuntimeError('gzip failed while exporting PLINK .raw')
            with gzip.open(packed, 'rb') as check:
                if not check.read(1): raise RuntimeError('Empty compressed PLINK output')
            packed.replace(output)
        finally:
            if keeper is not None: os.close(keeper)
            reader.close()
            for proc in (plink_proc, gzip_proc):
                if proc is not None and proc.poll() is None:
                    proc.terminate()
                    try: proc.wait(timeout=10)
                    except subprocess.TimeoutExpired: proc.kill(); proc.wait()
            log = Path(str(prefix)+'.log')
            if log.exists(): shutil.copyfile(log, str(output)+'.log')


def main():
    def stop(signum, frame): raise SystemExit(128+signum)
    for sig in (signal.SIGHUP, signal.SIGTERM): signal.signal(sig, stop)
    p = argparse.ArgumentParser()
    p.add_argument('--output', required=True)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    if not command or '--out' in command: p.error('supply PLINK command without --out')
    export(a.output, command)


if __name__ == '__main__': main()
