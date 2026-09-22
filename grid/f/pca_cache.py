#!/usr/bin/env python3
"""Check a projected-PC cache without loading the cohort into memory."""
import csv
import gzip
import sys
import math


def valid(path, n_pc):
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rt', newline='') as stream:
            reader = csv.reader(stream, delimiter='\t')
            header = next(reader)
            if not {'IID', '#IID', 'eid'}.intersection(header):
                return False
            if not {f'PC{i}' for i in range(1, n_pc + 1)}.issubset(header):
                return False
            count = 0
            idc=next(header.index(k) for k in ('IID','#IID','eid') if k in header)
            pcs=[header.index(f'PC{i}') for i in range(1,n_pc+1)]
            seen=set()
            for row in reader:
                if len(row) != len(header):
                    return False
                if row[idc] in ('', 'NA', 'NaN', 'nan') or row[idc] in seen or any(not math.isfinite(float(row[j])) for j in pcs):return False
                seen.add(row[idc])
                count += 1
            return count > 0
    except (OSError, EOFError, StopIteration, UnicodeError,ValueError):
        return False


if __name__ == '__main__':
    sys.exit(0 if valid(sys.argv[1], int(sys.argv[2])) else 1)
