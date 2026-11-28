"""Translate coordinate score IDs to target PVAR IDs without changing alleles/weights."""
import argparse
import csv
import subprocess
from collections import defaultdict
from contextlib import contextmanager


@contextmanager
def open_pvar(path):
    if str(path).endswith('.zst'):
        proc = subprocess.Popen(['zstd', '-dc', '--', str(path)], stdout=subprocess.PIPE, text=True)
        try:
            yield proc.stdout
        finally:
            proc.stdout.close()
            if proc.wait():
                raise RuntimeError('PVAR decompression failed')
    else:
        with open(path) as handle:
            yield handle


def key(chrom, pos, a, b):
    return (chrom.removeprefix('chr'), int(pos), tuple(sorted((a.upper(), b.upper()))))


def match_ids(score, pvar, output):
    with open(score) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        rows = list(reader)
    wanted = defaultdict(list)
    for i, row in enumerate(rows):
        parts = row[0].split(':')
        if len(parts) == 4 and parts[1].isdigit():
            wanted[key(*parts)].append(i)
    matches = defaultdict(set)
    if wanted:
        with open_pvar(pvar) as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 5:
                    continue
                chrom, pos, snp, ref, alt = fields[:5]
                if snp == '.' or ',' in alt:
                    continue
                for i in wanted.get(key(chrom, pos, ref, alt), ()):
                    if rows[i][1].upper() in (ref.upper(), alt.upper()):
                        matches[i].add(snp)
    coordinate_rows = {i for indices in wanted.values() for i in indices}
    retained = []
    for i, row in enumerate(rows):
        if i in coordinate_rows:
            if len(matches[i]) != 1:
                continue
            row[0] = next(iter(matches[i]))
        retained.append(row)
    ids = [row[0] for row in retained]
    if len(ids) != len(set(ids)):
        raise ValueError('Multiple score rows map to the same target ID; refusing duplicate weights')
    with open(output, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(retained)
    print(f'PGS ID matching: {len(rows)} input, {len(retained)} retained; '
          f'{len(coordinate_rows)} coordinate IDs checked against target PVAR', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('score')
    parser.add_argument('pvar')
    parser.add_argument('output')
    args = parser.parse_args()
    match_ids(args.score, args.pvar, args.output)
