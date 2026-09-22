#!/usr/bin/env python3
"""Store normalized GWAS beside the source; stage independent working copies."""
import argparse
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from prepare_sumstats import preparation_signature

FORMAT_VERSION = 1
COLUMNS = ['SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N', 'EAF', 'CHR', 'BP']


def cache_paths(source):
    source = Path(source).resolve()
    stem = source.name[:-3] if source.name.endswith('.gz') else source.name
    prefix = source.parent / (stem + '.csx.sumstats')
    return Path(str(prefix) + '.gz'), Path(str(prefix) + '.json')


def read_metadata(path):
    try:
        data = json.loads(Path(path).read_text())
        return data if isinstance(data, dict) else None
    except (OSError, ValueError):
        return None


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def validate_table(table, metadata):
    """Check the entire gzip stream and row count before adopting/publishing it."""
    with gzip.open(table, 'rt') as stream:
        if stream.readline().rstrip('\r\n').split('\t') != COLUMNS:
            raise ValueError(f'Not a normalized GWAS table: {table}')
        rows = sum(1 for line in stream if line.strip())
    if rows <= 0 or rows != metadata.get('kept_rows'):
        raise ValueError(f'Incomplete normalized GWAS: {table} ({rows} rows)')
    return sha256(table)


def cached_metadata(table, metadata, key):
    data = read_metadata(metadata)
    if (not data or data.get('preparation_signature') != key
            or data.get('cache_format') != FORMAT_VERSION):
        return None
    try:
        stat = table.stat()
        if (stat.st_size != data['table_size'] or stat.st_mtime_ns != data['table_mtime_ns']
                or stat.st_size == 0):
            return None
        with gzip.open(table, 'rt') as stream:
            if stream.readline().rstrip('\r\n').split('\t') != COLUMNS:
                return None
    except (OSError, EOFError, KeyError, ValueError):
        return None
    return data


def legacy_candidates(work, trait, pop, key, source):
    if work is None:
        return []
    found = []
    for path in sorted((Path(work) / trait).glob(f'*/sumstats/{pop}.json')):
        data = read_metadata(path)
        table = path.with_suffix('.tsv.gz')
        if (data and data.get('preparation_signature') == key
                and data.get('input') == str(source) and table.is_file()):
            found.append((table, path, data))
    return found


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream, indent=2)
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def stage_working_copy(table, data, output, metadata):
    """Keep per-run phi/N annotations and readers separate from the shared cache."""
    output, metadata = Path(output), Path(metadata)
    if output.resolve() == table.resolve():
        raise ValueError('Working output must differ from the permanent cache')
    output.parent.mkdir(parents=True, exist_ok=True)
    old = read_metadata(metadata)
    same = (old and old.get('preparation_signature') == data['preparation_signature']
            and old.get('table_sha256') == data['table_sha256'] and output.is_file()
            and output.stat().st_size == data['table_size']
            and output.stat().st_mtime_ns == data['table_mtime_ns'])
    if not same:
        fd, temporary = tempfile.mkstemp(prefix='.' + output.name + '.', dir=output.parent)
        os.close(fd)
        try:
            shutil.copy2(table, temporary)
            os.replace(temporary, output)
        finally:
            Path(temporary).unlink(missing_ok=True)
    run_data = dict(data, output=str(output.resolve()), permanent_output=str(table))
    atomic_json(metadata, run_data)


def ensure_prepared(source, snpinfo, trait, pop, *, work=None, output=None, metadata=None,
                    chunk=500000, replace=False, migrate_only=False, remove_legacy=False):
    source, snpinfo = Path(source).resolve(), Path(snpinfo).resolve()
    trait, pop = trait.lower(), pop.upper()
    if bool(output) != bool(metadata):
        raise ValueError('Working output and metadata must be supplied together')
    if remove_legacy and not migrate_only:
        raise ValueError('Removing old tables requires --migrate-only')
    table, info = cache_paths(source)
    key = preparation_signature(source, snpinfo, trait, pop)
    lock_path = info.with_suffix('.lock')
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        data = None if replace else cached_metadata(table, info, key)
        candidates = legacy_candidates(work, trait, pop, key, source)
        action = 'SKIP'
        if data is None:
            chosen = None
            if not replace:
                for old_table, old_info, old_data in candidates:
                    try:
                        digest = validate_table(old_table, old_data)
                    except (OSError, EOFError, ValueError):
                        continue
                    chosen = (old_table, old_data, digest)
                    break
            if chosen is None and migrate_only:
                print(f'MISSING preprocessing {trait}.{pop}: no complete matching table', flush=True)
                return None
            with tempfile.TemporaryDirectory(prefix='.' + table.name + '.', dir=table.parent) as temp:
                temporary = Path(temp) / 'sumstats.tsv.gz'
                temporary_info = Path(temp) / 'sumstats.json'
                if chosen:
                    old_table, data, digest = chosen
                    shutil.copy2(old_table, temporary)
                    if sha256(temporary) != digest:
                        raise ValueError(f'Cache copy verification failed: {old_table}')
                    action = 'MIGRATED'
                else:
                    print(f'RUN preprocessing {trait}.{pop}: {source}', flush=True)
                    subprocess.run([
                        sys.executable, str(Path(__file__).with_name('prepare_sumstats.py')),
                        '--input', str(source), '--output', str(temporary),
                        '--metadata', str(temporary_info), '--snpinfo', str(snpinfo),
                        '--trait', trait, '--pop', pop, '--chunk', str(chunk),
                    ], check=True)
                    data = read_metadata(temporary_info)
                    if not data or data.get('preparation_signature') != key:
                        raise ValueError('Preprocessing inputs changed during preparation')
                    digest = validate_table(temporary, data)
                    action = 'SAVED'
                if preparation_signature(source, snpinfo, trait, pop) != key:
                    raise ValueError('Preprocessing inputs changed during cache publication')
                data = dict(data)
                # These settings belong to a particular inference run, not normalization.
                for name in ('inference_phi', 'n_gwas_used', 'permanent_output'):
                    data.pop(name, None)
                stat = temporary.stat()
                data.update(output=str(table), cache_format=FORMAT_VERSION,
                            table_sha256=digest, table_size=stat.st_size,
                            table_mtime_ns=stat.st_mtime_ns)
                # The JSON is the commit marker; interrupted publication is never a cache hit.
                info.unlink(missing_ok=True)
                os.replace(temporary, table)
                atomic_json(info, data)
        if remove_legacy:
            for old_table, old_info, old_data in candidates:
                if old_table.resolve() == table.resolve():
                    continue
                if old_table.stat().st_size == data['table_size'] and sha256(old_table) == data['table_sha256']:
                    old_table.unlink()
                    # Retain per-run metadata (including phi/N) as provenance.
                    old_data['output'] = str(table)
                    old_data['permanent_output'] = str(table)
                    atomic_json(old_info, old_data)
        if output:
            stage_working_copy(table, data, output, metadata)
        print(f'{action} preprocessing {trait}.{pop}: {table}', flush=True)
        return table, info, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('input', 'snpinfo', 'trait', 'pop'):
        parser.add_argument('--' + key, required=True)
    for key in ('work', 'output', 'metadata'):
        parser.add_argument('--' + key)
    parser.add_argument('--chunk', type=int, default=500000)
    parser.add_argument('--replace', choices=('TRUE', 'FALSE'), default='FALSE')
    parser.add_argument('--migrate-only', action='store_true')
    parser.add_argument('--remove-legacy', action='store_true')
    args = parser.parse_args()
    if args.chunk < 1:
        parser.error('--chunk must be positive')
    ensure_prepared(args.input, args.snpinfo, args.trait, args.pop, work=args.work,
                    output=args.output, metadata=args.metadata, chunk=args.chunk,
                    replace=args.replace == 'TRUE', migrate_only=args.migrate_only,
                    remove_legacy=args.remove_legacy)


if __name__ == '__main__':
    main()
