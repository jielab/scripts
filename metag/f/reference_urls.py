#!/usr/bin/env python3
"""Export existing RefSeq selections as disjoint aria2 input files; no downloads."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from urllib.parse import quote, urlsplit
from reference_selection import activate, frozen_summary, load_profile, FOLDER

LIBRARIES = ('archaea', 'bacteria', 'viral', 'human', 'fungi', 'plant', 'protozoa',
             'vertebrate_mammalian', 'vertebrate_other', 'invertebrate', 'UniVec_Core')
PARTITIONS = {
    'ref_microbes_univec.url': ('archaea', 'bacteria', 'viral', 'fungi', 'protozoa', 'UniVec_Core'),
    'ref_plants.url': ('plant',),
    'ref_vertebrates.url': ('human', 'vertebrate_mammalian', 'vertebrate_other'),
    'ref_invertebrates.url': ('invertebrate',),
}


def selected(summary, library, db):
    digest = hashlib.sha256()
    header = None
    rows = {}
    with summary.open('rb', buffering=8 << 20) as stream:
        for raw in stream:
            digest.update(raw)
            line = raw.decode('utf-8').rstrip('\r\n')
            if line.startswith('#assembly_accession'):
                header = line[1:].split('\t')
            if not line or line.startswith('#'):
                continue
            if header is None:
                raise ValueError(f'Missing assembly summary header: {summary}')
            fields = line.split('\t')
            if len(fields) != len(header):
                raise ValueError(f'Incomplete assembly summary row: {summary}')
            row = dict(zip(header, fields))
            if row['assembly_level'].lower() not in ('chromosome', 'complete genome'):
                continue
            if row['ftp_path'] == 'na':
                continue
            remote = urlsplit(row['ftp_path'].rstrip('/'))
            if remote.hostname != 'ftp.ncbi.nlm.nih.gov' or remote.query or remote.fragment:
                raise ValueError(f'Unexpected NCBI assembly URL: {row["ftp_path"]}')
            path = PurePosixPath(remote.path.lstrip('/'))
            if '..' in path.parts or not str(path).startswith('genomes/all/GCF/'):
                raise ValueError(f'Unsafe/non-RefSeq path: {path}')
            path = path / (path.name + '_genomic.fna.gz')
            target = db / 'library' / library / str(path)
            if any(c in str(target) for c in '\r\n\t'):
                raise ValueError(f'Unsafe aria2 option path: {target!s}')
            rows[str(path)] = dict(library=library, target=str(target),
                                  url='https://ftp.ncbi.nlm.nih.gov/' + quote(str(path), safe='/%'),
                                  bases=int(row['genome_size']))
    if not rows:
        raise ValueError(f'No chromosome/complete-genome references in {summary}')
    return rows, digest.hexdigest()


def export(args):
    db, out = args.db.resolve(), args.out_dir.resolve()
    summaries = (args.summary_dir or db.parent / 'reference_size_plan').resolve()
    if getattr(args, 'microbes_simple', False):
        activate(db, summaries)
    profile = load_profile(db)
    partitions = dict(PARTITIONS)
    if profile:
        partitions = {('ref_microbes.url' if name == 'ref_microbes_univec.url' else name): libraries
                      for name, libraries in partitions.items()}
    if any((db / name).exists() for name in ('hash.k2d', 'opts.k2d', 'taxo.k2d')):
        raise ValueError('An index already exists; do not download into a frozen/partial index.')
    records, sources, skipped, per_library = [], {}, [], {}
    for library in LIBRARIES:
        lib = db / 'library' / library
        marker = db / 'metag_build_state' / (library + '.complete')
        if marker.is_file() and all((lib / name).is_file() and (lib / name).stat().st_size > 0
                                   for name in ('library.fna', 'prelim_map.txt')):
            skipped.append(library)
            continue
        if library == 'UniVec_Core':
            items = [dict(library=library, target=str(lib / library),
                          url='https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core', bases=1000000)]
        else:
            selected_snapshot = frozen_summary(db, library)
            summary = selected_snapshot or lib / 'assembly_summary.txt'
            if not summary.is_file():
                summary = summaries / (library + '.assembly_summary.txt')
            if not summary.is_file():
                raise ValueError(f'Missing local assembly summary: {summary}')
            rows, digest = selected(summary, library, db)
            manifest = lib / 'manifest.txt'
            if manifest.is_file():
                paths = [line.strip() for line in manifest.read_text().splitlines() if line.strip()]
                if len(set(paths)) != len(paths) or set(paths) != set(rows):
                    original = db / FOLDER / 'originals' / f'{library}.manifest.txt'
                    if selected_snapshot is None or not original.is_file() or manifest.read_bytes() != original.read_bytes():
                        raise ValueError(f'Existing manifest does not match summary selection: {manifest}')
            sources[library] = {'path': str(summary), 'sha256': digest}
            items = list(rows.values())
        records.extend(items)
        per_library[library] = len(items)
    targets = [r['target'] for r in records]
    if len(targets) != len(set(targets)):
        raise ValueError('Overlapping download destinations')
    if not records:
        raise ValueError('All source libraries are already complete; no URL lists needed.')
    # Biological categories make each queue's contents explicit.
    queues = [[r for r in sorted(records, key=lambda r: r['target']) if r['library'] in libraries]
              for libraries in partitions.values()]
    if sum(map(len, queues)) != len(records):
        raise ValueError('Reference categories are not covered exactly once by the four partitions.')
    outputs = {}
    for filename, rows in zip(partitions, queues):
        lines = ['# Generated from existing metag RefSeq selections; aria2 input format.\n']
        for row in rows:
            target = Path(row['target'])
            lines.append(f'{row["url"]}\n  dir={target.parent}\n  out={target.name}\n')
        outputs[out / filename] = ''.join(lines)
    existing = set(out.glob('ref_*.url'))
    if existing - outputs.keys():
        raise ValueError('Old partitions are present; use a new output directory.')
    for path, text in outputs.items():
        if path.exists() and path.read_text() != text:
            raise ValueError(f'Refusing to overwrite a different download list: {path}')
    out.mkdir(parents=True, exist_ok=True)
    for path, text in outputs.items():
        if not path.exists():
            with path.open('x') as stream:
                stream.write(text)
    receipt = dict(created_utc=datetime.now(timezone.utc).isoformat(), db=str(db),
                   selection='RefSeq chromosome + complete_genome, same as build-db',
                   profile=profile['profile'] if profile else 'all_eligible_assemblies',
                   completed_libraries_skipped=skipped, counts=per_library, total_files=len(records),
                   partitioning='biological categories; file counts and sizes are not balanced',
                   sources=sources, integrity='HTTPS prefetch only; build-db still verifies NCBI MD5 before use.',
                   parts=[dict(file=str(path), count=len(rows), estimated_bases=sum(r['bases'] for r in rows),
                               libraries=list(partitions[path.name]),
                               sha256=hashlib.sha256(text.encode()).hexdigest())
                          for (path, text), rows in zip(outputs.items(), queues)])
    receipt_path = out / 'metag_reference_urls.json'
    temp = receipt_path.with_suffix('.json.partial')
    temp.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
    temp.replace(receipt_path)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--db', type=Path, default=Path('/mnt/f/metag/ref_kraken2_ictv/kraken2_broad_v1'))
    p.add_argument('--out-dir', type=Path, default=Path('/mnt/f/metag/ref_kraken2_ictv'))
    p.add_argument('--summary-dir', type=Path, help='Existing supplementary assembly summary snapshots')
    p.add_argument('--microbes-simple', action='store_true',
                   help='Freeze one eligible assembly per microbial species; only before affected libraries are processed')
    p.add_argument('--parts', type=int, choices=(4,), default=4,
                   help='Fixed four biological groups: microbes/UniVec, plants, vertebrates, invertebrates')
    args = p.parse_args()
    try:
        export(args)
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
