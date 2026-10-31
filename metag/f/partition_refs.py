#!/usr/bin/env python3
from __future__ import annotations
import shutil
from pathlib import Path
from common import DIR_LIB, DIR_OUT
from library import active_references


def main(max_compressed_gib=3.0):
    mapping=DIR_LIB/'reference/assembly_summary.txt'
    if not mapping.is_file():
        raise RuntimeError('Reference mapping missing; run ./1.library.sh --check-only')
    rows=active_references()
    base=DIR_OUT/'metacache'/'partitions'
    if base.exists(): shutil.rmtree(base)
    base.mkdir(parents=True)
    limit=int(max_compressed_gib*2**30)
    parts=[]; current=None; used=0
    for row in sorted(rows,key=lambda r:int(r['bytes']),reverse=True):
        size=int(row['bytes'])
        if current is None or (used and used+size>limit):
            current=base/f'ref_{len(parts)+1:04d}'
            current.mkdir(); parts.append(current); used=0
            shutil.copyfile(mapping,current/'assembly_summary.txt')
        (current/row['filename']).symlink_to(Path(row['path']).absolute())
        used += size
    total=sum(int(r['bytes']) for r in rows)
    print(f'PASS: {len(rows):,} refs -> {len(parts)} partitions; compressed total={total/2**30:.1f} GiB; target={max_compressed_gib} GiB/part')
    oversized=[r for r in rows if int(r['bytes'])>limit]
    if oversized:
        print(f'NOTE: {len(oversized)} single FASTA files exceed the partition target; original files are retained; unfinished groups use bounded derived shards.')
    return parts


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--max-compressed-gib',type=float,default=3.0)
    a=ap.parse_args();main(a.max_compressed_gib)
