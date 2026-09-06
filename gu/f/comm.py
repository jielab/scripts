#!/usr/bin/env python3
"""Compatibility helpers for GU TSV files with very wide fields."""
from __future__ import annotations

import csv
import gzip
import os
import sys
import tempfile
from pathlib import Path


def resolve_tsv_path(path: Path) -> Path:
    """Prefer compressed unfiltered inventories, with legacy plain TSV fallback."""
    path = Path(path)
    if path.name.endswith('unfiltered.tsv.gz'):
        plain, compressed = path.with_suffix(''), path
    elif path.name.endswith('unfiltered.tsv'):
        plain, compressed = path, path.with_name(path.name + '.gz')
    else:
        return path
    if compressed.is_file():
        return compressed
    return plain if plain.is_file() else path


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    path = resolve_tsv_path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv_rows(path: Path, fieldnames, rows) -> None:
    """Stream TSV output atomically; unfiltered inventories are always gzip."""
    path = Path(path)
    if path.name.endswith('unfiltered.tsv'):
        path = path.with_name(path.name + '.gz')
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f'.{path.name}.', suffix='.tmp', dir=path.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(tmp, 'wt', encoding='utf-8', newline='') as handle:
            fields = list(fieldnames)
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t',
                                    lineterminator='\n', extrasaction='ignore')
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, '') for key in fields})
        tmp.replace(path)
        if path.name.endswith('unfiltered.tsv.gz'):
            # Retire a previous plain output only after the new gzip is closed
            # and atomically published. Failed writes preserve both old forms.
            path.with_suffix('').unlink(missing_ok=True)
    finally:
        tmp.unlink(missing_ok=True)


def enable_wide_csv_fields() -> int:
    """Raise the process-wide CSV field limit to the platform maximum."""
    limit = sys.maxsize
    while limit > 0:
        try:
            csv.field_size_limit(limit)
            return limit
        except OverflowError:
            # Some Python builds expose a C long smaller than sys.maxsize.
            limit //= 10
    raise RuntimeError("unable to configure a usable CSV field-size limit")


CHROM_LENGTHS = {
    "37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560,
    },
    "38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895,
    },
}

"""Expand normalized BED loci while retaining an auditable core/analysis map."""
import argparse
import re
from pathlib import Path

def parse_bp(value:str)->int:
    m=re.fullmatch(r'\s*([0-9]+)\s*(bp|k|kb|m|mb)?\s*',value,re.I)
    if not m:raise argparse.ArgumentTypeError('use a non-negative size such as 0, 100kb, or 1mb')
    scale={None:1,'bp':1,'k':1000,'kb':1000,'m':1000000,'mb':1000000}[m.group(2).lower() if m.group(2) else None]
    return int(m.group(1))*scale

def expand_loci_main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True,type=Path); ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--map',required=True,type=Path); ap.add_argument('--flank',default='100kb',type=parse_bp)
    ap.add_argument('--build',required=True,choices=['37','38']); args=ap.parse_args()
    rows=[]
    for line_no,line in enumerate(args.input.read_text().splitlines(),1):
        if not line.strip() or line.lstrip().startswith('#'):continue
        fields=line.split('\t');
        if len(fields)<4:raise SystemExit(f'ERROR: normalized BED line {line_no} has fewer than four columns')
        chrom,start,end,locus=fields[:4]; start=int(start); end=int(end); chrom_len=CHROM_LENGTHS[args.build][chrom]
        analysis_start=max(0,start-args.flank); analysis_end=min(chrom_len,end+args.flank)
        rows.append((chrom,start,end,analysis_start,analysis_end,locus,args.flank))
    if not rows:raise SystemExit('ERROR: no loci to expand')
    args.output.parent.mkdir(parents=True,exist_ok=True); args.map.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(''.join(f'{c}\t{a0}\t{a1}\t{locus}\n' for c,_,_,a0,a1,locus,_ in rows))
    args.map.write_text('chr\tcore_start\tcore_end\tanalysis_start\tanalysis_end\tlocus_id\tflank_bp\n'+
                        ''.join(f'{c}\t{s}\t{e}\t{a0}\t{a1}\t{locus}\t{flank}\n' for c,s,e,a0,a1,locus,flank in rows))
    print(args.flank)


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print('Usage: comm.py expand-loci [options]')
    elif sys.argv.pop(1) == 'expand-loci':
        expand_loci_main()
    else:
        raise SystemExit('Unknown shared utility command')
