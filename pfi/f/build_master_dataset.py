#!/usr/bin/env python3


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import argparse
import glob
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import duckdb
import pyarrow.parquet as pq

PUBMED_COLUMNS = [
    'iso_journal', 'issn', 'publication_types', 'mesh_terms', 'keywords',
    'authors', 'first_author', 'last_author', 'author_count', 'affiliations',
]

PMC_COLUMNS = [
    'article_type', 'license', 'ref_count', 'table_count', 'fig_count',
    'intro_text', 'methods_text', 'results_text', 'discussion_text',
    'data_availability_text', 'ethics_text', 'trial_registration_text',
    'full_text',
]

def log(message: str) -> None:
    now = datetime.now().isoformat(timespec='seconds')
    print(f'[{now}] {message}', flush=True)

def sql_string(value: os.PathLike | str) -> str:
    return str(value).replace("'", "''")

def parquet_glob(path: Path) -> str:
    if path.is_dir() or path.suffix != '.parquet':
        return sql_string(path / '**' / '*.parquet')
    return sql_string(path)

def parquet_files_exist(path: Path) -> bool:
    if path.is_file() and path.suffix == '.parquet':
        return True
    return path.exists() and any(path.rglob('*.parquet'))

def parquet_row_count(path: Path) -> int:
    """Read row counts from Parquet footers without scanning file contents."""
    if path.is_file() and path.suffix == '.parquet':
        return int(pq.ParquetFile(path).metadata.num_rows)
    if not path.exists():
        return 0
    total = 0
    for file_path in path.rglob('*.parquet'):
        total += int(pq.ParquetFile(file_path).metadata.num_rows)
    return total

def success_marker(path: Path) -> Path:
    return path / '_SUCCESS'

def row_count_marker(path: Path) -> Path:
    return path / '_ROW_COUNT'

def read_row_count_marker(path: Path) -> int | None:
    marker = row_count_marker(path)
    if not marker.exists():
        return None
    try:
        return int(marker.read_text(encoding='utf-8').strip())
    except ValueError:
        return None

def mark_success(path: Path, rows: int | None = None) -> None:
    if path.is_dir():
        success_marker(path).write_text(datetime.now().isoformat(timespec='seconds') + '\n')
        if rows is not None:
            row_count_marker(path).write_text(f'{int(rows)}\n', encoding='utf-8')

def source_dir_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for child in path.iterdir() if child.is_dir() and child.name.startswith('source_'))

def reusable_source_dataset(path: Path, expected_sources: int | None = None) -> bool:
    if success_marker(path).exists():
        return True
    return path.exists() and expected_sources is not None and source_dir_count(path) >= expected_sources

def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()

def connect(temp_dir: Path | None, threads: int, memory_limit: str | None = None) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(f'PRAGMA threads={int(threads)}')
    con.execute('SET preserve_insertion_order=false')
    if memory_limit:
        con.execute(f"SET memory_limit = '{sql_string(memory_limit)}'")
    if temp_dir:
        temp_dir.mkdir(parents=True, exist_ok=True)
        con.execute(f"SET temp_directory = '{sql_string(temp_dir)}'")
    return con

def copy_partitioned(con: duckdb.DuckDBPyConnection, query: str, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"""
        COPY ({query})
        TO '{sql_string(out_dir)}'
        (FORMAT PARQUET, COMPRESSION ZSTD, PARTITION_BY (bucket_partition),
         OVERWRITE_OR_IGNORE true)
    """)

def normalize_pubmed(args, pubmed_dir: Path) -> None:
    if args.reuse_pubmed and parquet_files_exist(pubmed_dir):
        log(f'Reusing existing PubMed hash-bucket partitions: {pubmed_dir}')
        n = read_row_count_marker(pubmed_dir)
        if n is None:
            n = parquet_row_count(pubmed_dir)
        log(f'pubmed deduplicated rows: {n:,}')
        return
    log('Building deduplicated PubMed hash-bucket partitions...')
    con = connect(Path(args.temp_dir) if args.temp_dir else None, args.threads, args.memory_limit)
    pubmed_glob = sql_string(args.pubmed_glob)
    query = f"""
        SELECT *,
               cast(hash(article_key) % {int(args.buckets)} AS INTEGER) AS bucket_partition
        FROM (
            SELECT * EXCLUDE(rn),
                   coalesce(nullif(pmid, ''), nullif(pmcid, ''), nullif(doi, ''), nullif(title, '')) AS article_key
            FROM (
                SELECT *, row_number() OVER (
                    PARTITION BY coalesce(nullif(pmid, ''), nullif(pmcid, ''), nullif(doi, ''), nullif(title, ''))
                    ORDER BY year DESC NULLS LAST
                ) AS rn
                FROM read_parquet('{pubmed_glob}', union_by_name=true)
            )
            WHERE rn = 1
        )
    """
    copy_partitioned(con, query, pubmed_dir)
    n = parquet_row_count(pubmed_dir)
    log(f'pubmed deduplicated rows: {n:,}')
    mark_success(pubmed_dir, n)
    con.close()

def list_parquet_files(pattern: str) -> list[str]:
    files = glob.glob(pattern, recursive=True)
    return sorted(f for f in files if f.endswith('.parquet'))

def normalize_one_pmc_file(task: tuple[int, str, str, str | None, int, int, str | None]) -> tuple[int, int]:
    index, path, out_root, temp_dir, threads, buckets, memory_limit = task
    source_dir = Path(out_root) / f'source_{index:05d}'
    remove_path(source_dir)
    con = connect(Path(temp_dir) if temp_dir else None, threads, memory_limit)
    pmc_path = sql_string(path)
    query = f"""
        SELECT *,
               cast(hash(article_key) % {int(buckets)} AS INTEGER) AS bucket_partition
        FROM (
            SELECT *,
                   coalesce(nullif(pmid, ''), nullif(pmcid, ''), nullif(doi, ''), nullif(title, '')) AS article_key
            FROM read_parquet('{pmc_path}', union_by_name=true)
        )
    """
    copy_partitioned(con, query, source_dir)
    n = parquet_row_count(source_dir)
    con.close()
    return index, n

def normalize_pmc(args, pmc_dir: Path) -> None:
    files = list_parquet_files(args.pmc_glob)
    if not files:
        raise FileNotFoundError(f'No PMC parquet files matched: {args.pmc_glob}')
    log(f'PMC source parquet files: {len(files):,}')
    if args.reuse_pmc and reusable_source_dataset(pmc_dir, expected_sources=len(files)):
        log(f'Reusing existing PMC OA hash-bucket partitions: {pmc_dir}')
        n = read_row_count_marker(pmc_dir)
        if n is not None:
            log(f'pmc rows written: {n:,}')
        return
    log('Building PMC OA hash-bucket partitions by source Parquet file...')
    remove_path(pmc_dir)
    pmc_dir.mkdir(parents=True, exist_ok=True)
    tasks = [
        (i, path, str(pmc_dir), args.temp_dir, args.threads, args.buckets, args.memory_limit)
        for i, path in enumerate(files)
    ]
    total = 0
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=max(1, int(args.jobs))) as executor:
        futures = [executor.submit(normalize_one_pmc_file, task) for task in tasks]
        done = 0
        for future in as_completed(futures):
            index, n = future.result()
            total += n
            done += 1
            if done <= 10 or done % 100 == 0 or done == len(tasks):
                elapsed = max(time.monotonic() - started, 1.0)
                files_per_min = done / elapsed * 60.0
                remaining = max(len(tasks) - done, 0)
                eta_min = remaining / files_per_min if files_per_min > 0 else 0.0
                log(
                    f'PMC partitioned files: {done:,}/{len(tasks):,}; '
                    f'last source={index:05d}; rows total so far={total:,}; '
                    f'rate={files_per_min:.1f} files/min; eta={eta_min:.1f} min'
                )
    mark_success(pmc_dir, total)
    log(f'pmc rows written: {total:,}')

def bucket_sets(pubmed_dir: Path, pmc_dir: Path, bucket_count: int) -> tuple[list[int], set[int]]:
    pubmed_buckets = set()
    if pubmed_dir.exists():
        for child in pubmed_dir.iterdir():
            if child.is_dir() and child.name.startswith('bucket_partition='):
                raw = child.name.split('=', 1)[1]
                try:
                    pubmed_buckets.add(int(raw))
                except ValueError:
                    pass

    # PMC is written as source_*/bucket_partition=*/part files. Scanning every
    # source directory can take more than an hour on Windows-mounted drives, so
    # rely on the configured hash-bucket count once the dataset exists.
    pmc_buckets = set(range(int(bucket_count))) if parquet_files_exist(pmc_dir) else set()
    buckets = sorted(pubmed_buckets | pmc_buckets)
    if not buckets and (pubmed_dir.exists() or pmc_dir.exists()):
        buckets = list(range(int(bucket_count)))
    return buckets, pmc_buckets

def null_select(alias: str, columns: list[str]) -> list[str]:
    return [f'NULL AS {col}' for col in columns]

def select_master(p_alias: str | None, m_alias: str | None) -> str:
    def coalesce_col(col: str) -> str:
        if p_alias and m_alias:
            return f'coalesce({p_alias}.{col}, {m_alias}.{col}) AS {col}'
        if p_alias:
            return f'{p_alias}.{col} AS {col}'
        return f'{m_alias}.{col} AS {col}'

    parts = [
        coalesce_col('pmid'),
        coalesce_col('pmcid'),
        coalesce_col('doi'),
        coalesce_col('year'),
        coalesce_col('title'),
        coalesce_col('abstract'),
        coalesce_col('journal'),
    ]
    if p_alias:
        parts.extend(f'{p_alias}.{col} AS {col}' for col in PUBMED_COLUMNS)
    else:
        parts.extend(null_select('p', PUBMED_COLUMNS))
    if m_alias:
        parts.extend(f'{m_alias}.{col} AS {col}' for col in PMC_COLUMNS)
        parts.append(
            f"CASE WHEN coalesce({m_alias}.full_text, {m_alias}.intro_text, {m_alias}.methods_text, "
            f"{m_alias}.results_text, {m_alias}.discussion_text, {m_alias}.data_availability_text, "
            f"{m_alias}.ethics_text, {m_alias}.trial_registration_text) IS NOT NULL "
            "THEN 1 ELSE 0 END AS has_pmc_full_text"
        )
    else:
        parts.extend(null_select('m', PMC_COLUMNS))
        parts.append('0 AS has_pmc_full_text')
    return ',\n            '.join(parts)

def join_one_bucket(task: tuple[int, str, str, str, str | None, int, bool]) -> tuple[int, int]:
    bucket, pubmed_dir, pmc_dir, out_dir, temp_dir, threads, m_exists = task
    p_dir = Path(pubmed_dir) / f'bucket_partition={bucket}'
    b_out = Path(out_dir) / f'bucket_partition={bucket}'
    remove_path(b_out)
    b_out.mkdir(parents=True, exist_ok=True)
    con = connect(Path(temp_dir) if temp_dir else None, threads)

    p_exists = p_dir.exists()
    m_glob = sql_string(Path(pmc_dir) / '*' / f'bucket_partition={bucket}' / '*.parquet')
    if p_exists and m_exists:
        query = f"""
            SELECT
            {select_master('p', 'm')}
            FROM read_parquet('{parquet_glob(p_dir)}', union_by_name=true) p
            FULL OUTER JOIN read_parquet('{m_glob}', union_by_name=true) m
            USING(article_key)
        """
    elif p_exists:
        query = f"""
            SELECT
            {select_master('p', None)}
            FROM read_parquet('{parquet_glob(p_dir)}', union_by_name=true) p
        """
    elif m_exists:
        query = f"""
            SELECT
            {select_master(None, 'm')}
            FROM read_parquet('{m_glob}', union_by_name=true) m
        """
    else:
        con.close()
        return bucket, 0

    con.execute(f"""
        COPY ({query})
        TO '{sql_string(b_out / 'part.parquet')}'
        (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    n = parquet_row_count(b_out / 'part.parquet')
    con.close()
    return bucket, n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pubmed_glob', required=True)
    ap.add_argument('--pmc_glob', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--work_dir', default='/mnt/d/analysis/pfi/work/parquet/06_build')
    ap.add_argument('--temp_dir', default=None)
    ap.add_argument('--jobs', type=int, default=2)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--buckets', type=int, default=128)
    ap.add_argument('--memory_limit', default=None)
    ap.add_argument('--progress_every_buckets', type=int, default=8)
    ap.add_argument('--fresh_work', action='store_true')
    ap.add_argument('--reuse_pubmed', action='store_true')
    ap.add_argument('--reuse_pmc', action='store_true')

    # Backward-compatible no-op arguments accepted by older shell scripts.
    ap.add_argument('--db', default=None)
    ap.add_argument('--fresh_db', action='store_true')
    ap.add_argument('--delete_db_after', action='store_true')
    args = ap.parse_args()

    work_dir = Path(args.work_dir)
    pubmed_dir = work_dir / 'pubmed_by_bucket'
    old_pubmed_dir = work_dir / 'pubmed_by_year'
    if old_pubmed_dir.exists() and not pubmed_dir.exists():
        old_pubmed_dir.rename(pubmed_dir)
    pmc_dir = work_dir / 'pmc_by_bucket'
    out_dir = Path(args.out)
    Path(args.temp_dir).mkdir(parents=True, exist_ok=True) if args.temp_dir else None
    work_dir.mkdir(parents=True, exist_ok=True)

    if args.fresh_work:
        if not args.reuse_pubmed:
            remove_path(pubmed_dir)
        if not args.reuse_pmc:
            remove_path(pmc_dir)
        remove_path(out_dir)

    normalize_pubmed(args, pubmed_dir)
    normalize_pmc(args, pmc_dir)

    buckets, pmc_buckets = bucket_sets(pubmed_dir, pmc_dir, args.buckets)
    log(f'Joining PubMed and PMC OA by hash buckets: {len(buckets)} bucket groups; jobs={args.jobs}; duckdb_threads_per_job={args.threads}')
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = [
        (bucket, str(pubmed_dir), str(pmc_dir), str(out_dir), args.temp_dir, args.threads, bucket in pmc_buckets)
        for bucket in buckets
    ]
    total = 0
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=max(1, int(args.jobs))) as executor:
        futures = [executor.submit(join_one_bucket, task) for task in tasks]
        done = 0
        progress_every = max(1, int(args.progress_every_buckets))
        for future in as_completed(futures):
            bucket, n = future.result()
            total += n
            done += 1
            if done <= 2 or done % progress_every == 0 or done == len(tasks):
                avg_rows = total / done
                est_total = int(round(avg_rows * len(tasks)))
                elapsed = max(time.monotonic() - started, 1.0)
                buckets_per_min = done / elapsed * 60.0
                remaining = max(len(tasks) - done, 0)
                eta_min = remaining / buckets_per_min if buckets_per_min > 0 else 0.0
                log(
                    f'joined buckets: {done:,}/{len(tasks):,}; '
                    f'last bucket_partition={bucket}; last rows={n:,}; '
                    f'rows so far={total:,}; estimated final rows={est_total:,}; '
                    f'rate={buckets_per_min:.1f} buckets/min; eta={eta_min:.1f} min'
                )

    log(f'articles_master rows written: {total:,}')
    log(f'articles_master dataset: {out_dir}')
    mark_success(out_dir, total)

if __name__ == '__main__':
    main()
