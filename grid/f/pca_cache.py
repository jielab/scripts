#!/usr/bin/env python3
"""Check a projected-PC cache without loading the cohort into memory."""
import csv
import gzip
import sys


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
            for row in reader:
                if len(row) != len(header):
                    return False
                count += 1
            return count > 0
    except (OSError, EOFError, StopIteration, UnicodeError):
        return False


if __name__ == '__main__':
    sys.exit(0 if valid(sys.argv[1], int(sys.argv[2])) else 1)
