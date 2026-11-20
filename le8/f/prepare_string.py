#!/usr/bin/env python3
"""Convert already downloaded human STRING physical links to an offline table.

No network access. Supply matching-version human protein.info and
protein.physical.links files (plain text or gzip); do not use functional links.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile


def text_file(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def convert(info, links, directory):
    names = {}
    with text_file(info) as f:
        reader = csv.DictReader(f, delimiter='\t')
        id_col = next((k for k in ('#string_protein_id', 'string_protein_id', 'protein_external_id') if k in (reader.fieldnames or [])), None)
        if not id_col or 'preferred_name' not in (reader.fieldnames or []):
            raise ValueError('protein.info requires string_protein_id and preferred_name columns')
        for row in reader:
            if row[id_col].startswith('9606.') and row['preferred_name']:
                names[row[id_col]] = row['preferred_name']
    if not names:
        raise ValueError('No human (9606) identifiers in protein.info')
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / 'string_physical.tsv'
    if target.exists():
        raise FileExistsError(f'{target} exists; use a new snapshot directory to replace a database')
    count = 0
    fd, temporary = tempfile.mkstemp(dir=directory)
    try:
        with os.fdopen(fd, 'w') as out, text_file(links) as src:
            writer = csv.writer(out, delimiter='\t', lineterminator='\n')
            writer.writerow(['preferredName_A', 'preferredName_B', 'score'])
            header = src.readline().strip().split()
            required = ['protein1', 'protein2', 'combined_score']
            if not set(required).issubset(header):
                raise ValueError('physical.links requires protein1 protein2 combined_score columns')
            indexes = [header.index(k) for k in required]
            for line in src:
                if not line.strip(): continue
                fields = line.split()
                a, b, score = [fields[i] for i in indexes]
                score = float(score)
                if not 0 <= score <= 1000:
                    raise ValueError('Expected STRING combined_score on a 0..1000 scale')
                if score >= 400 and a in names and b in names and names[a].upper() != names[b].upper():
                    writer.writerow([names[a], names[b], format(score / 1000, '.6g')])
                    count += 1
        if not count:
            raise ValueError('No eligible human physical edges; check matching database versions')
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)
    h = hashlib.sha256()
    with target.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    (directory / 'string_source.json').write_text(json.dumps(dict(
        info=str(Path(info).resolve()), physical_links=str(Path(links).resolve()),
        file=target.name, sha256=h.hexdigest(), minimum_score=400, human_edges=count), indent=2) + '\n')
    return count


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--info', type=Path, required=True)
    p.add_argument('--physical-links', type=Path, required=True)
    p.add_argument('--cache', type=Path, required=True)
    args = p.parse_args()
    print('Offline physical edges:', convert(args.info, args.physical_links, args.cache))
