#!/usr/bin/env python3
"""Remove only managed reference FASTA files after a complete, current audit."""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import errno
import fcntl
import json
import os
from pathlib import Path
import stat
import sys
import time
import uuid

from reference_layout import DATA_DIRS
from taxonomy import atomic_json, sha256, snapshot


def plain_path(path):
    """Do not follow symlinks, including linked source/library directories."""
    for part in (path, *path.parents):
        if part.is_symlink():
            raise RuntimeError(f'Reference cleanup refuses symlink: {part}')


def managed(db):
    return (db / 'metag_build_state/config.txt').is_file()


def validate(db, mapping):
    state = db / 'metag_build_state'
    for name in ('all.complete', 'build.complete', 'audit.complete'):
        if not (state / name).is_file():
            raise RuntimeError(f'Reference build is not complete: {state / name}')
    config = state / 'config.txt'
    entries = [line.partition('=')[2].split() for line in config.read_text().splitlines()
               if line.startswith('libraries=')]
    if (len(entries) != 1 or not entries[0] or len(set(entries[0])) != len(entries[0])
            or any(lib not in DATA_DIRS for lib in entries[0])):
        raise RuntimeError('Cleanup requires a known, unambiguous library configuration')
    libraries = entries[0]
    for lib in libraries:
        if not (state / f'{lib}.complete').is_file():
            raise RuntimeError(f'Reference library is not complete: {lib}')
    for name in ('hash.k2d', 'opts.k2d', 'taxo.k2d'):
        path = db / name
        plain_path(path)
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f'Missing/empty index: {path}')
    audit = json.loads((db / 'metag.audit.json').read_text())
    if (audit.get('schema') != 'v2' or audit.get('missing_groups') != []
            or audit.get('missing_virus_classes') != []
            or audit.get('virus_types_tested') is not True
            or audit.get('exploratory_virus_coverage', False)):
        raise RuntimeError('Cleanup requires successful base-group AND DNA/RNA virus coverage audits')
    current = snapshot(db)
    if current != audit.get('snapshot'):
        raise RuntimeError('Index/taxonomy changed since audit; refusing cleanup')
    plain_path(mapping)
    mapping_hash = sha256(mapping)
    if mapping_hash != (audit.get('virus_mapping') or {}).get('mapping_sha256'):
        raise RuntimeError('Virus mapping changed since audit; refusing cleanup')
    inspect = db / 'metag.inspect.tsv'
    if not inspect.is_file() or sha256(inspect) != audit.get('inspect_sha256'):
        raise RuntimeError('Missing/changed index inspection evidence; run a fresh db-audit')
    return libraries, dict(snapshot=current, mapping_sha256=mapping_hash,
                           config_sha256=sha256(config), audit_sha256=sha256(db / 'metag.audit.json'))


def file_identity(path):
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise RuntimeError(f'Cleanup requires a regular file: {path}')
    return dict(bytes=info.st_size, mtime_ns=info.st_mtime_ns,
                inode=info.st_ino, device=info.st_dev)


def candidates(db, libraries):
    """Explicit generated filenames only; never match FASTQ or arbitrary *.gz."""
    files, directories = {}, set()

    def add(path, kind):
        plain_path(path)
        if not path.exists():
            return
        for suffix in ('.aria2', '.tmp'):
            if path.with_name(path.name + suffix).exists():
                raise RuntimeError(f'Unfinished reference download/processing: {path}')
        files[str(path)] = dict(path=str(path), kind=kind, **file_identity(path))

    def walk(folder):
        plain_path(folder)
        if not folder.exists():
            return
        with os.scandir(folder) as scan:
            for entry in scan:
                path = Path(entry.path)
                if entry.is_symlink():
                    raise RuntimeError(f'Reference cleanup refuses symlink: {path}')
                if entry.is_dir(follow_symlinks=False):
                    directories.add(path)
                    walk(path)
                elif entry.name.endswith(('.aria2', '.tmp')):
                    raise RuntimeError(f'Unfinished reference download/processing: {path}')
                elif entry.name.endswith('_genomic.fna'):
                    add(path, 'processed_genome')
                elif entry.name.endswith('_genomic.fna.gz'):
                    add(path, 'reference_download')

    for lib in libraries:
        folder = db / 'library' / lib
        add(folder / 'library.fna', 'merged_library')
        genomes = folder / 'genomes'
        walk(genomes)
        if genomes.exists():
            directories.add(genomes)
        # The shared layout stores downloads outside the database. Only inspect
        # each known leaf directory, never the raw/split trees or all of annot/.
        if db.parent.name == '0.metadata':
            source = db.parent.parent / DATA_DIRS[lib]
            plain_path(source)
            if source.exists():
                with os.scandir(source) as scan:
                    for entry in scan:
                        if entry.name.endswith(('.aria2', '.tmp')):
                            raise RuntimeError(f'Unfinished reference download: {entry.path}')
                        if entry.name.endswith('_genomic.fna.gz'):
                            add(Path(entry.path), 'reference_download')
    return sorted(files.values(), key=lambda item: item['path']), directories


def remove_file(item):
    path = Path(item['path'])
    plain_path(path)
    if file_identity(path) != {k: item[k] for k in ('bytes', 'mtime_ns', 'inode', 'device')}:
        raise RuntimeError(f'Reference file changed after cleanup planning: {path}')
    path.unlink()


def cleanup(db, mapping, execute=False):
    db, mapping = Path(os.path.abspath(db)), Path(os.path.abspath(mapping))
    plain_path(db)
    if not managed(db):
        # Imported databases may have unrelated files; do not assume ownership.
        print('[metag] SKIP reference cleanup: database is not managed by metag', flush=True)
        return dict(status='unmanaged', deleted_files=0, deleted_bytes=0)
    with ExitStack() as stack:
        for name in ('.metag_build.lock', '.metag.audit.lock'):
            path = db / name
            plain_path(path)
            handle = stack.enter_context(path.open('a'))
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError('Reference build/audit/cleanup is still running; no cleanup performed') from error
        libraries, evidence = validate(db, mapping)
        print('[metag] START reference cleanup: checking generated FASTA and reference .fna.gz files', flush=True)
        files, directories = candidates(db, libraries)
        total = sum(item['bytes'] for item in files)
        print(f'[metag] Reference cleanup plan: {len(files):,} files, {total / 2**30:.2f} GiB', flush=True)
        if not execute:
            return dict(status='planned', files=files, total_bytes=total)
        if not files:
            print('[metag] SKIP reference cleanup: no reference FASTA files remain', flush=True)
            return dict(status='empty', deleted_files=0, deleted_bytes=0)
        # Recheck evidence after potentially slow traversal, before any deletion.
        if validate(db, mapping)[1] != evidence:
            raise RuntimeError('Reference evidence changed during cleanup planning')
        log_dir = db / 'metag_cleanup'
        plain_path(log_dir)
        log_dir.mkdir(exist_ok=True)
        token = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '.' + uuid.uuid4().hex[:8]
        plan = log_dir / (token + '.plan.jsonl')
        events = log_dir / (token + '.events.jsonl')
        summary = dict(schema='metag-reference-cleanup-v1', status='running',
                       started_utc=datetime.now(timezone.utc).isoformat(), evidence=evidence,
                       plan=str(plan), events=str(events), planned_files=len(files),
                       planned_bytes=total, deleted_files=0, deleted_bytes=0,
                       policy='Generated reference FASTA and reference _genomic.fna.gz only; sample FASTQ retained')
        with plan.open('x') as stream:
            for item in files:
                stream.write(json.dumps(item, ensure_ascii=False) + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        report = log_dir / (token + '.summary.json')
        atomic_json(report, summary)
        try:
            with events.open('x', buffering=1) as stream:
                last = time.monotonic()
                for item in files:
                    remove_file(item)
                    summary['deleted_files'] += 1
                    summary['deleted_bytes'] += item['bytes']
                    stream.write(json.dumps(dict(path=item['path'], bytes=item['bytes'], status='deleted'),
                                            ensure_ascii=False) + '\n')
                    if time.monotonic() - last >= 60:
                        print(f'[metag] Reference cleanup: {summary["deleted_files"]:,}/{len(files):,} files, '
                              f'{summary["deleted_bytes"] / 2**30:.2f} GiB removed', flush=True)
                        stream.flush()
                        os.fsync(stream.fileno())
                        last = time.monotonic()
                stream.flush()
                os.fsync(stream.fileno())
            # Prune only empty generated genome directories, never metadata trees.
            for path in sorted(directories, key=lambda p: len(p.parts), reverse=True):
                plain_path(path)
                try:
                    path.rmdir()
                except OSError as error:
                    if error.errno not in (errno.ENOTEMPTY, errno.ENOENT, errno.EEXIST):
                        raise
            if validate(db, mapping)[1] != evidence:
                raise RuntimeError('Reference evidence changed during cleanup')
            summary['status'] = 'complete'
        except BaseException as error:
            summary.update(status='interrupted', error=str(error))
            raise
        finally:
            summary['finished_utc'] = datetime.now(timezone.utc).isoformat()
            atomic_json(report, summary)
        atomic_json(db / 'metag.cleanup.json', summary)
        print(f'[metag] DONE reference cleanup: removed {summary["deleted_files"]:,} files, '
              f'{summary["deleted_bytes"] / 2**30:.2f} GiB; receipt={report}', flush=True)
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', '--dir-db', type=Path, required=True)
    parser.add_argument('--virus-map', type=Path, required=True)
    parser.add_argument('--execute', action='store_true', help='Delete after validation; default only lists planned files')
    args = parser.parse_args()
    try:
        result = cleanup(args.db, args.virus_map, args.execute)
        if not args.execute:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except (OSError, RuntimeError, ValueError) as error:
        print(f'[metag] Reference cleanup stopped: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
