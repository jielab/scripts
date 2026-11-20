#!/usr/bin/env python3
"""Run one independent command phase; a failed phase stops the downstream script."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
from workflow_cpu import worker_cores


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('commands', type=Path)
    p.add_argument('--workers', type=int, required=True)
    a=p.parse_args()
    commands=json.loads(a.commands.read_text())
    if a.workers < 1: p.error('--workers must be positive')
    env={**os.environ, 'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    def run(command):
        subprocess.run(['bash','-euo','pipefail','-c',command],env=env,check=True)
    with ThreadPoolExecutor(max_workers=min(a.workers,worker_cores())) as pool:
        for _ in pool.map(run,commands): pass

if __name__ == '__main__': main()
