#!/usr/bin/env python3
"""Publish canonical MAHA diet RDS files; preserve separate DOSE caches."""
import argparse
import os
from pathlib import Path
import shutil
import tempfile


def publish(source, destination):
    source, destination = Path(source), Path(destination)
    files = sorted(source.glob('diet*.rds'))
    if not {'diet0.rds', 'diet.rds'} <= {p.name for p in files}:
        raise FileNotFoundError('MAHA scoring must complete diet0.rds and diet.rds before publication')
    if source.resolve() == destination.resolve():
        raise ValueError('Cache and publication directory must differ')
    destination.mkdir(parents=True, exist_ok=True)
    staged = []
    try:
        for path in files:
            fd, name = tempfile.mkstemp(prefix='.diet-publish-', dir=destination)
            os.close(fd)
            tmp = Path(name)
            staged.append((tmp, destination / path.name))
            shutil.copy2(path, tmp)
        for tmp, target in staged:
            os.replace(tmp, target)
    finally:
        for tmp, _ in staged:
            tmp.unlink(missing_ok=True)
    print(f'Published {len(files)} MAHA diet RDS files to {destination}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    publish(args.source, args.destination)
