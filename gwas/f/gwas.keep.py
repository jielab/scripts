#!/usr/bin/env python3
"""Restrict REGENIE's allocated sample rows to the prepared complete cases."""
import argparse
import os
from pathlib import Path
import tempfile


def write_keep(phenotype, output, keep=None):
    allowed = None
    if keep:
        allowed = set()
        with open(keep) as stream:
            for line in stream:
                row = line.split()
                if not row or [x.lstrip('#') for x in row] == ['FID', 'IID']:
                    continue
                if len(row) != 2:
                    raise ValueError('REGENIE --keep requires two columns: FID IID')
                allowed.add(tuple(row))
    selected, seen = [], set()
    with open(phenotype) as stream:
        header = stream.readline().split()
        fi, ii = header.index('FID'), header.index('IID')
        for line in stream:
            row = line.split()
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError('Malformed prepared phenotype row')
            key = row[fi], row[ii]
            if key in seen:
                raise ValueError('Duplicate prepared phenotype sample')
            seen.add(key)
            if allowed is None or key in allowed:
                selected.append(key)
    if not selected:
        raise ValueError('No complete-case samples remain after --keep intersection')
    output = Path(output)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=output.parent, delete=False) as stream:
            temporary = stream.name
            for fid, iid in selected:
                stream.write(f'{fid}\t{iid}\n')
        os.replace(temporary, output)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
    return len(selected)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phe', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--keep')
    args = parser.parse_args()
    print(f'REGENIE complete-case keep: {write_keep(args.phe, args.out, args.keep)} samples')
