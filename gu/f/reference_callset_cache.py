#!/usr/bin/env python3
"""Build and validate GU's persistent normalized published-callset database."""
from __future__ import annotations

import argparse
import os
import shutil
import sqlite3
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from normalize_results import (
    REFERENCE_COLS,
    REFERENCE_RAW_COLS,
    collapse_reference_callsets,
    iter_reference_callsets,
)

CACHE_FORMAT='1'


def callset_files(source):
    return sorted(
        path for path in source.rglob('*')
        if path.is_file() and path.stat().st_size > 0
        and (path.name.lower().endswith(('.bed','.tsv','.bed.gz','.tsv.gz')))
    )


def require_callset_files(source):
    files=callset_files(source)
    if not files:raise ValueError(f'no BED/TSV callsets found under {source}')
    return files


def inspect_cache(cache,quiet=False,full=True):
    try:
        if not cache.is_file():raise ValueError(f'cache is missing: {cache}')
        con=sqlite3.connect(f'file:{cache}?mode=ro',uri=True)
        try:
            meta=dict(con.execute('SELECT key,value FROM metadata'))
            if meta.get('cache_format')!=CACHE_FORMAT:raise ValueError('cache format changed')
            raw=con.execute('SELECT COUNT(*) FROM reference_callsets_raw').fetchone()[0]
            union=con.execute('SELECT COUNT(*) FROM reference_callsets').fetchone()[0]
            if full:
                check=con.execute('PRAGMA quick_check').fetchone()
                if not check or check[0]!='ok':raise ValueError(f'SQLite quick_check failed: {check}')
            if raw<1 or union<1:raise ValueError(f'cache tables are empty: raw={raw} union={union}')
        finally:con.close()
        if not quiet:print(f'Reference callset database OK: raw={raw} union={union} database={cache}')
        return True
    except (OSError,sqlite3.Error,ValueError) as exc:
        if not quiet:print(f'ERROR: missing or invalid reference callset database: {exc}',file=os.sys.stderr)
        return False


def build_cache(source,cache,dataset_id):
    files=require_callset_files(source)
    cache.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp_name=tempfile.mkstemp(prefix='gu-reference-callsets-',suffix='.sqlite'); os.close(fd)
    tmp=Path(tmp_name); staged=cache.with_name(f'.{cache.name}.part.{os.getpid()}')
    con=None; raw=0
    try:
        con=sqlite3.connect(tmp)
        con.executescript('''
          PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL; PRAGMA temp_store=FILE;
          CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
          CREATE TABLE reference_callsets_raw(
            dataset_id TEXT,population TEXT,genome_build TEXT,chr TEXT,start INTEGER,end INTEGER,
            source_class TEXT,reference_role TEXT,raw_file TEXT,reference_sample_id TEXT,reference_super_population TEXT
          );
          CREATE TABLE reference_callsets(
            dataset_id TEXT,population TEXT,genome_build TEXT,chr TEXT,start INTEGER,end INTEGER,
            source_class TEXT,reference_role TEXT,raw_file TEXT
          );
        ''')
        insert_sql='INSERT INTO reference_callsets_raw VALUES(?,?,?,?,?,?,?,?,?,?,?)'
        for frame in iter_reference_callsets(source,dataset_id):
            con.executemany(insert_sql,frame[REFERENCE_RAW_COLS].itertuples(index=False,name=None)); raw+=len(frame)
        con.commit(); collapse_reference_callsets(con)
        con.executescript('''
          CREATE INDEX idx_reference_raw_sample_region ON reference_callsets_raw(reference_sample_id,genome_build,source_class,chr,start,end);
          CREATE INDEX idx_reference_callsets_region ON reference_callsets(genome_build,source_class,chr,start,end);
          CREATE INDEX idx_reference_callsets_population ON reference_callsets(dataset_id,population,genome_build);
          ANALYZE;
        ''')
        union=con.execute('SELECT COUNT(*) FROM reference_callsets').fetchone()[0]
        metadata={
            'cache_format':CACHE_FORMAT,
            'dataset_id':dataset_id,'source_dir':str(source.resolve()),'source_files':str(len(files)),
            'raw_rows':str(raw),'union_rows':str(union),'created_utc':datetime.now(timezone.utc).isoformat(),
        }
        con.executemany('INSERT INTO metadata(key,value) VALUES(?,?)',metadata.items()); con.commit()
        check=con.execute('PRAGMA quick_check').fetchone()
        if not check or check[0]!='ok':raise RuntimeError(f'SQLite quick_check failed: {check}')
        con.close(); con=None
        if staged.exists():staged.unlink()
        shutil.copy2(tmp,staged); os.replace(staged,cache)
    finally:
        if con is not None:con.close()
        if tmp.exists():tmp.unlink()
        if staged.exists():staged.unlink()
    print(f'Reference callset database built: files={len(files)} raw={raw} union={union} database={cache}')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('action',choices=('prepare','check'))
    ap.add_argument('--source',required=True,type=Path)
    ap.add_argument('--cache',required=True,type=Path)
    ap.add_argument('--dataset-id',default='AS3_1KG')
    ap.add_argument('--fast',action='store_true',help='Skip the full SQLite page scan; schema and row counts are still checked')
    args=ap.parse_args()
    if args.action=='check':raise SystemExit(0 if inspect_cache(args.cache,full=not args.fast) else 2)
    if inspect_cache(args.cache,quiet=True,full=False):
        inspect_cache(args.cache); return
    if not args.source.is_dir():raise SystemExit(f'ERROR: published callset directory is missing: {args.source}')
    build_cache(args.source,args.cache,args.dataset_id)
    if not inspect_cache(args.cache):raise SystemExit(2)


if __name__=='__main__':main()
