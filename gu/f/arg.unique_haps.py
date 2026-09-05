#!/usr/bin/env python3
"""Drop all records at repeated positions from sorted PLINK Oxford HAPS.

Keep only one pending row in memory. Do not select an arbitrary allele from
split multiallelic sites. Replace the input only after a successful full scan.
"""
import argparse
import json
import os
from pathlib import Path


def clean(path, chromosome):
    path = Path(path)
    temporary = path.with_name(path.name + '.unique.next')
    total = kept = duplicate_sites = duplicate_records = 0
    previous = None
    pending = None
    count = 0
    try:
        with path.open() as source, temporary.open('w', newline='\n') as output:
            for number, line in enumerate(source, 1):
                if not line.strip():
                    continue
                fields = line.split(maxsplit=5)
                if len(fields) != 6 or fields[0].removeprefix('chr') != str(chromosome):
                    raise ValueError(f'HAPS line {number}: malformed row or wrong chromosome')
                position = int(fields[2])
                if position < 1 or (previous is not None and position < previous):
                    raise ValueError(f'HAPS line {number}: invalid/unsorted position {position} after {previous}')
                total += 1
                if position == previous:
                    if count == 1:
                        duplicate_sites += 1
                        duplicate_records += 1
                    duplicate_records += 1
                    count += 1
                    pending = None
                else:
                    if pending is not None:
                        output.write(pending)
                        kept += 1
                    previous, pending, count = position, line, 1
            if pending is not None:
                output.write(pending)
                kept += 1
        if kept < 2:
            raise ValueError(f'Too few unique-position HAPS variants: {kept}')
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return dict(input_variants=total, retained_variants=kept,
                duplicate_positions=duplicate_sites, removed_variants=duplicate_records,
                policy='exclude-all-repeated-positions')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--haps', required=True)
    parser.add_argument('--chr', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    result = clean(args.haps, args.chr)
    Path(args.out).write_text(json.dumps(result, indent=2) + '\n')
    print('HAPS position QC: ' + json.dumps(result))


if __name__ == '__main__':
    main()
