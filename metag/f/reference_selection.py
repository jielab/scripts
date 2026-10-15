#!/usr/bin/env python3
"""Freeze a species-deduplicated microbial reference selection before library processing."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import uuid
from urllib.parse import urlsplit

PROFILE = 'microbial_simple_v1'
LIBRARIES = ('bacteria', 'fungi', 'protozoa')
FOLDER = 'metag_reference_selection'


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            result.update(block)
    return result.hexdigest()


def load_profile(db):
    path = Path(db) / FOLDER / 'profile.json'
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    if data.get('profile') != PROFILE or set(data.get('libraries', {})) != set(LIBRARIES):
        raise ValueError('Unrecognized frozen reference selection')
    config = Path(db) / 'metag_build_state/config.txt'
    if config.exists() and f'selection_sha256={digest(path)}' not in config.read_text().splitlines():
        raise ValueError('Frozen selection does not match the build configuration')
    return data


def frozen_summary(db, library):
    data = load_profile(db)
    if data is None or library not in LIBRARIES:
        return None
    path = Path(db) / FOLDER / f'{library}.assembly_summary.txt'
    if digest(path) != data['libraries'][library]['selected_sha256']:
        raise ValueError(f'Frozen assembly selection changed: {path}')
    return path


def summary_for_k2(args):
    """Called before k2 changes directories or overwrites any library metadata."""
    path = frozen_summary(Path(args.db).resolve(), args.library)
    if path is not None and (args.protein or args.assembly_source != 'refseq'
                            or set(args.assembly_levels) != {'chromosome', 'complete_genome'}):
        raise ValueError('Frozen microbial selection requires nucleotide RefSeq chromosome/complete_genome settings')
    return str(path) if path else None


def choose_species(summary):
    best, header, candidates, total_bases = {}, None, 0, 0
    source_hash = hashlib.sha256()
    with summary.open('rb', buffering=8 << 20) as stream:
        for raw in stream:
            source_hash.update(raw)
            line = raw.decode('utf-8').rstrip('\r\n')
            if line.startswith('#assembly_accession'):
                header = line[1:].split('\t')
            if not line or line.startswith('#'):
                continue
            if header is None:
                raise ValueError(f'Missing assembly summary header: {summary}')
            values = line.split('\t')
            if len(values) != len(header):
                raise ValueError(f'Incomplete assembly row: {summary}')
            row = dict(zip(header, values))
            if row['assembly_level'].lower() not in ('complete genome', 'chromosome') or row['ftp_path'] == 'na':
                continue
            if not row['species_taxid'].isdigit() or int(row['species_taxid']) <= 0:
                raise ValueError(f'Missing species taxid for {row["assembly_accession"]}')
            remote = urlsplit(row['ftp_path'].rstrip('/'))
            if remote.hostname != 'ftp.ncbi.nlm.nih.gov' or not remote.path.startswith('/genomes/all/GCF/'):
                raise ValueError('Expected an NCBI RefSeq assembly URL')
            candidates += 1
            total_bases += int(row['genome_size'])
            rank = (row['refseq_category'] == 'reference genome',
                    row['refseq_category'] == 'representative genome',
                    row['assembly_level'] == 'Complete Genome', row['seq_rel_date'], row['assembly_accession'])
            species = row['species_taxid']
            if species not in best or rank > best[species][0]:
                best[species] = (rank, row, line)
    if not best:
        raise ValueError(f'No eligible assemblies: {summary}')
    chosen = sorted(best.values(), key=lambda x: x[1]['assembly_accession'])
    text = '#' + '\t'.join(header) + '\n' + ''.join(item[2] + '\n' for item in chosen)
    info = dict(source_sha256=source_hash.hexdigest(), candidate_assemblies=candidates,
                selected_assemblies=len(chosen), candidate_species=len(best),
                selected_species=len(best), candidate_bases=total_bases,
                selected_bases=sum(int(item[1]['genome_size']) for item in chosen),
                selected_sha256=hashlib.sha256(text.encode()).hexdigest())
    return text, info


def signature(db):
    data = load_profile(db)
    if data is None:
        return ''
    for library in LIBRARIES:
        frozen_summary(db, library)
    return f'selection_profile={PROFILE}\nselection_sha256={digest(Path(db) / FOLDER / "profile.json")}'


def activate(db, summaries=None):
    db = Path(db).resolve()
    summaries = Path(summaries).resolve() if summaries else db.parent / 'reference_size_plan'
    db.mkdir(parents=True, exist_ok=True)
    with (db / '.metag_build.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if load_profile(db) is not None:
            signature(db)
            return
        if any((db / name).exists() for name in ('hash.k2d', 'opts.k2d', 'taxo.k2d')):
            raise ValueError('An index already exists; choose a NEW database for a different selection')
        for library in LIBRARIES:
            paths = [db / 'metag_build_state' / f'{library}.complete']
            paths += [db / 'library' / library / n for n in ('library.fna', 'prelim_map.txt')]
            if any(p.exists() for p in paths):
                raise ValueError(f'{library} has already been processed; choose a NEW database')
        config = db / 'metag_build_state/config.txt'
        previous_config = config.read_text().rstrip('\n') if config.exists() else None
        if previous_config and any(line.startswith('selection_') for line in previous_config.splitlines()):
            raise ValueError('Existing selection configuration is incomplete; inspect before retrying')
        staged = db / (FOLDER + '.partial.' + uuid.uuid4().hex)
        staged.mkdir()
        originals = staged / 'originals'
        originals.mkdir()
        info = dict(profile=PROFILE, created_utc=datetime.now(timezone.utc).isoformat(),
                    policy='One assembly per species_taxid within the existing complete-genome/chromosome candidates; '
                           'prefer NCBI reference, representative, complete genome, newest release, then accession.',
                    retained_completed_libraries=['archaea', 'viral'],
                    other_categories='Unchanged', libraries={})
        try:
            for library in LIBRARIES:
                source = db / 'library' / library / 'assembly_summary.txt'
                if not source.exists():
                    source = summaries / f'{library}.assembly_summary.txt'
                text, details = choose_species(source)
                (staged / f'{library}.assembly_summary.txt').write_text(text)
                original = originals / f'{library}.assembly_summary.txt'
                shutil.copyfile(source, original)
                if digest(original) != details['source_sha256']:
                    raise ValueError('Source summary changed while selecting references')
                details['original_summary'] = f'originals/{original.name}'
                info['libraries'][library] = details
                manifest = db / 'library' / library / 'manifest.txt'
                if manifest.exists():
                    shutil.copyfile(manifest, originals / f'{library}.manifest.txt')
            if previous_config is not None:
                (originals / 'build_config.txt').write_text(previous_config + '\n')
            profile = staged / 'profile.json'
            profile.write_text(json.dumps(info, ensure_ascii=False, indent=2) + '\n')
            checksum = digest(profile)
            staged.rename(db / FOLDER)
            if previous_config is not None:
                tmp = config.with_suffix('.selection.tmp')
                tmp.write_text(previous_config + f'\nselection_profile={PROFILE}\nselection_sha256={checksum}\n')
                tmp.replace(config)
        except BaseException:
            # Keep a published profile for diagnosis; only remove our unpublished metadata staging directory.
            if staged.exists():
                shutil.rmtree(staged)
            raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('activate', 'signature'))
    p.add_argument('--db', type=Path, required=True)
    p.add_argument('--summary-dir', type=Path)
    args = p.parse_args()
    try:
        if args.action == 'activate':
            activate(args.db, args.summary_dir)
        else:
            print(signature(args.db), end='')
    except (ValueError, OSError, KeyError) as error:
        p.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
