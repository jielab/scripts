"""Validate the shared male-only, one-haplotype-per-person ARG X contract."""
import argparse
import csv
from pathlib import Path


def pfile_prefix(directory, chrom):
    return Path(directory) / ('chrX.male' if chrom == 'X' else f'chr{chrom}')


def male_sample_ids(psam):
    with Path(psam).open() as handle:
        header = handle.readline().lstrip('#').split()
        if not {'IID', 'SEX'} <= set(header):
            raise ValueError(f'chrX.male PSAM requires IID and SEX: {psam}')
        rows = [dict(zip(header, line.split())) for line in handle if line.strip()]
    names = [row['IID'] for row in rows]
    if not names or len(set(names)) != len(names) or any(row.get('SEX') != '1' for row in rows):
        raise ValueError(f'chrX.male requires unique male samples (SEX=1): {psam}')
    return set(names)


def validate_x_map(sample_map, psam):
    males = male_sample_ids(psam)
    with Path(sample_map).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    names = [row['sample'] for row in rows]
    if (not names or not set(names) <= males or len(set(names)) != len(names)
            or any(row['haplotype'] != '1' for row in rows)):
        raise ValueError('chrX ARG must contain only chrX.male samples, each with haplotype=1; '
                         'rebuild chrX with arg.sh build --chr X --format trace')
    return len(names)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--psam', type=Path, required=True)
    parser.add_argument('--sample-map', type=Path)
    args = parser.parse_args()
    count = (validate_x_map(args.sample_map, args.psam) if args.sample_map
             else len(male_sample_ids(args.psam)))
    print(f'ARG chrX male-only sample check PASS: samples={count} source={args.psam}')
