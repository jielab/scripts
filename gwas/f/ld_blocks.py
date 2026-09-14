#!/usr/bin/env python3
"""Download the published GRCh37 LDetect blocks and record their provenance."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
from urllib.request import urlopen


def read_bed(text, omitted=None):
    rows = []
    for line in text.splitlines():
        fields = line.split()
        if not fields or fields[0].lower() in ('chr', 'chrom', '#chrom'):
            continue
        chromosome = fields[0].removeprefix('chr')
        if chromosome not in {str(i) for i in range(1, 23)}:
            raise ValueError('Unexpected LDetect chromosome: ' + chromosome)
        if 'None' in fields[1:3] and omitted is not None:
            omitted.append(line)
            continue
        start, end = map(int, fields[1:3])
        if not 0 <= start < end:
            raise ValueError('Invalid BED interval: ' + line)
        rows.append((int(chromosome), start, end))
    rows.sort()
    if {row[0] for row in rows} != set(range(1, 23)):
        raise ValueError('Expected all 22 autosomes')
    for left, right in zip(rows, rows[1:]):
        if left[0] == right[0] and left[2] > right[1]:
            raise ValueError('Overlapping BED intervals')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default='/mnt/i/refLD/block')
    args = parser.parse_args()
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    records = []
    for race in ('AFR', 'ASN', 'EUR'):
        url = f'https://api.bitbucket.org/2.0/repositories/nygcresearch/ldetect-data/src/master/{race}/fourier_ls-all.bed'
        with urlopen(url, timeout=60) as response:
            raw = response.read()
        omitted = []
        rows = read_bed(raw.decode(), omitted)
        # Upstream AFR chr11 contains two records with an unknown shared edge.
        # Omit these invalid intervals; do not invent a replacement breakpoint.
        if omitted:
            (root / f'{race}.37.omitted.tsv').write_text('\n'.join(omitted) + '\n')
        data = ''.join(f'chr{chrom}\t{start}\t{end}\n' for chrom, start, end in rows).encode()
        targets = ('EAS', 'SAS') if race == 'ASN' else (race,)
        # Do not overwrite a different locally curated block set on reruns.
        for target in targets:
            path = root / f'{target}.37.bed'
            if path.exists() and path.read_bytes() != data:
                raise FileExistsError('Existing BED differs: ' + str(path))
        stage = root / f'{race}.37.bed'
        if stage.exists() and stage.read_bytes() != data:
            raise FileExistsError('Existing BED differs: ' + str(stage))
        stage.write_bytes(data)
        if race == 'ASN':
            for target in targets:
                shutil.copyfile(stage, root / f'{target}.37.bed')
            stage.unlink()
        for target in targets:
            records.append([f'{target}.37.bed', '37', target, race, len(rows),
                            '0-based half-open; chromosomes 1-22', url,
                            hashlib.sha256(raw).hexdigest(), hashlib.sha256(data).hexdigest()])
            print(f'{target}.37.bed: {len(rows)} blocks; source={race}', flush=True)
    for race in ('AFR', 'EAS', 'EUR', 'SAS'):
        path = root / f'{race}.38.bed'
        if path.exists():
            rows = read_bed(path.read_text())
            records.append([path.name, '38', race, race, len(rows),
                            '0-based half-open; chromosomes 1-22',
                            'https://github.com/jmacdon/LDblocks_GRCh38',
                            '', hashlib.sha256(path.read_bytes()).hexdigest()])
    with (root / 'block_sources.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['file', 'grch', 'race', 'source_population', 'blocks', 'coordinates',
                         'source_url', 'download_sha256', 'file_sha256'])
        writer.writerows(records)
    (root / 'README.bplot.md').write_text(
        '# LD block references for bplot\n\n'
        'GRCh37: published LDetect fourier_ls-all.bed files, with BED coordinates retained.\n'
        'AFR: two upstream chr11 intervals contain a None endpoint and are omitted, '
        'recorded in AFR.37.omitted.tsv. No replacement boundary was inferred.\n'
        'EAS.37.bed and SAS.37.bed are identical copies of ASN, as requested; '
        'they are not independently inferred EAS/SAS block sets.\n\n'
        'GRCh38: existing pyrho files renamed to [race].38.bed without changing contents.\n'
        'Both sets contain autosomes 1-22 only. No HIS, ALL or chrX boundary set is supplied.\n'
        'See block_sources.tsv for source URLs, population mappings and checksums.\n')


if __name__ == '__main__':
    main()
