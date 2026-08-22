#!/usr/bin/env python3

import argparse
import os
import shutil
from collections import defaultdict
from pathlib import Path
from datetime import datetime

import duckdb
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


TEXT_SECTIONS = [
    ('title', 'TITLE'),
    ('abstract', 'ABSTRACT'),
    ('intro_text', 'INTRODUCTION'),
    ('methods_text', 'METHODS'),
    ('results_text', 'RESULTS'),
    ('discussion_text', 'DISCUSSION'),
    ('data_availability_text', 'DATA AVAILABILITY'),
    ('ethics_text', 'ETHICS'),
    ('trial_registration_text', 'TRIAL REGISTRATION'),
]

FULL_TEXT_MARKER_COLUMNS = [
    'full_text',
    'intro_text',
    'methods_text',
    'results_text',
    'discussion_text',
    'data_availability_text',
    'ethics_text',
    'trial_registration_text',
]

OUTPUT_META_COLUMNS = [
    'pmid',
    'pmcid',
    'doi',
    'year',
    'publication_types',
    'mesh_terms',
    'keywords',
]


def log(message):
    now = datetime.now().isoformat(timespec='seconds')
    print(f'[{now}] {message}', flush=True)


def sql_string(value):
    return str(value).replace("'", "''")


def parquet_source(path):
    p = Path(path)
    if p.is_dir() or p.suffix != '.parquet':
        return sql_string(p / '**' / '*.parquet')
    return sql_string(p)


def is_temp_parquet_path(file_path):
    return any(part.startswith('.year=') and '.tmp-' in part for part in file_path.parts)


def parquet_row_count(path, strict=True):
    p = Path(path)
    if not p.exists():
        return 0
    try:
        if p.is_file() and p.suffix == '.parquet':
            return int(pq.ParquetFile(p).metadata.num_rows)
        total = 0
        for file_path in p.rglob('*.parquet'):
            if is_temp_parquet_path(file_path):
                continue
            total += int(pq.ParquetFile(file_path).metadata.num_rows)
        return total
    except Exception:
        if strict:
            raise
        return -1


def blinded_select(full_text_expr, include_year=True):
    year_col = 'year,' if include_year else ''
    return f'''
            SELECT
                pmid, pmcid, doi, {year_col}
                publication_types, mesh_terms, keywords,
                CASE WHEN {full_text_expr} THEN 1 ELSE 0 END AS has_pmc_full_text,
                trim(concat_ws('\n\n',
                    CASE WHEN title IS NOT NULL THEN 'TITLE: ' || title ELSE NULL END,
                    CASE WHEN abstract IS NOT NULL THEN 'ABSTRACT: ' || abstract ELSE NULL END,
                    CASE WHEN intro_text IS NOT NULL THEN 'INTRODUCTION: ' || intro_text ELSE NULL END,
                    CASE WHEN methods_text IS NOT NULL THEN 'METHODS: ' || methods_text ELSE NULL END,
                    CASE WHEN results_text IS NOT NULL THEN 'RESULTS: ' || results_text ELSE NULL END,
                    CASE WHEN discussion_text IS NOT NULL THEN 'DISCUSSION: ' || discussion_text ELSE NULL END,
                    CASE WHEN data_availability_text IS NOT NULL THEN 'DATA AVAILABILITY: ' || data_availability_text ELSE NULL END,
                    CASE WHEN ethics_text IS NOT NULL THEN 'ETHICS: ' || ethics_text ELSE NULL END,
                    CASE WHEN trial_registration_text IS NOT NULL THEN 'TRIAL REGISTRATION: ' || trial_registration_text ELSE NULL END
                )) AS blinded_text
            FROM m
    '''


def hive_year_value(value):
    if value is None:
        return '__HIVE_DEFAULT_PARTITION__'
    try:
        number = int(value)
        if float(value) == float(number):
            return str(number)
    except (TypeError, ValueError):
        pass
    return str(value).replace('/', '_').replace('\\', '_')


def year_filter(value):
    if value is None:
        return 'year IS NULL'
    return f'year = {int(value)}'


def full_text_expr_sql():
    return (
        'coalesce(full_text, intro_text, methods_text, results_text, discussion_text, '
        'data_availability_text, ethics_text, trial_registration_text) IS NOT NULL'
    )


def where_clause_sql(full_text_only):
    if full_text_only:
        return full_text_expr_sql()
    return 'title IS NOT NULL OR abstract IS NOT NULL OR methods_text IS NOT NULL OR results_text IS NOT NULL'


def connect_duckdb(args):
    con = duckdb.connect()
    con.execute(f'PRAGMA threads={int(args.threads)}')
    con.execute('SET preserve_insertion_order = false')
    if args.memory_limit:
        con.execute(f"SET memory_limit = '{sql_string(args.memory_limit)}'")
    if args.temp_dir:
        Path(args.temp_dir).mkdir(parents=True, exist_ok=True)
        con.execute(f"SET temp_directory = '{sql_string(args.temp_dir)}'")
    return con


def duckdb_partitions(args, where_clause):
    con = connect_duckdb(args)
    try:
        con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet('{parquet_source(args.master)}', union_by_name=true)")
        return con.execute(f'''
            SELECT year, count(*) AS n
            FROM m
            WHERE {where_clause}
            GROUP BY year
            ORDER BY year NULLS LAST
        ''').fetchall()
    finally:
        con.close()


def clean_stale_temp_dirs(out_path):
    if not out_path.exists() or out_path.is_file():
        return
    for child in out_path.iterdir():
        if child.is_dir() and child.name.startswith('.year=') and '.tmp-' in child.name:
            shutil.rmtree(child)


def partition_dir(out_path, year):
    return out_path / f'year={hive_year_value(year)}'


def existing_partition_status(out_path, partitions):
    status = {}
    for year, expected_rows in partitions:
        part_dir = partition_dir(out_path, year)
        actual_rows = parquet_row_count(part_dir, strict=False)
        status[year] = (actual_rows, int(expected_rows))
    return status


def partition_years_to_write(out_path, partitions, resume):
    if not resume:
        return {year for year, _ in partitions}
    status = existing_partition_status(out_path, partitions)
    years = set()
    for year, (actual_rows, expected_rows) in status.items():
        if actual_rows == expected_rows and expected_rows > 0:
            log(f'Skipping complete partition year={year}; rows={actual_rows:,}')
        else:
            if actual_rows > 0:
                log(
                    f'Rewriting incomplete/stale partition year={year}; '
                    f'existing_rows={actual_rows:,}; expected_rows={expected_rows:,}'
                )
            years.add(year)
    return years


def open_pyarrow_dataset(path):
    p = Path(path)
    if p.is_dir():
        return ds.dataset(p, format='parquet', partitioning='hive')
    return ds.dataset(p, format='parquet')


def field_type(schema, name, default=pa.string()):
    try:
        return schema.field(name).type
    except KeyError:
        return default


def output_schema(input_schema, include_year):
    fields = []
    for name in OUTPUT_META_COLUMNS:
        if name == 'year' and not include_year:
            continue
        fields.append(pa.field(name, field_type(input_schema, name)))
    fields.append(pa.field('has_pmc_full_text', pa.int64()))
    fields.append(pa.field('blinded_text', pa.string()))
    return pa.schema(fields)


def valid_any_expr(column_names):
    expr = None
    for name in column_names:
        item = ds.field(name).is_valid()
        expr = item if expr is None else (expr | item)
    return expr


def pyarrow_filter_expr(schema_names, full_text_only, years_to_write):
    if full_text_only:
        text_columns = [c for c in FULL_TEXT_MARKER_COLUMNS if c in schema_names]
    else:
        text_columns = [c for c in ['title', 'abstract', 'methods_text', 'results_text'] if c in schema_names]
    expr = valid_any_expr(text_columns) if text_columns else None

    if years_to_write is not None and 'year' in schema_names:
        non_null_years = [int(y) for y in years_to_write if y is not None]
        year_expr = None
        if non_null_years:
            year_expr = ds.field('year').isin(non_null_years)
        if any(y is None for y in years_to_write):
            null_expr = ds.field('year').is_null()
            year_expr = null_expr if year_expr is None else (year_expr | null_expr)
        if year_expr is not None:
            expr = year_expr if expr is None else (expr & year_expr)
    return expr


def values_for_column(cols, name, n):
    return cols[name] if name in cols else [None] * n


def row_has_any_value(cols, names, index):
    for name in names:
        if name in cols and cols[name][index] is not None:
            return True
    return False


def build_blinded_text(cols, index):
    parts = []
    for name, label in TEXT_SECTIONS:
        if name not in cols:
            continue
        value = cols[name][index]
        if value is not None:
            parts.append(f'{label}: {value}')
    return '\n\n'.join(parts).strip()


def batch_to_tables_by_year(batch, schema, full_text_only, years_to_write):
    cols = batch.to_pydict()
    n = batch.num_rows
    if n == 0:
        return {}

    if full_text_only:
        required_text_columns = [c for c in FULL_TEXT_MARKER_COLUMNS if c in cols]
    else:
        required_text_columns = [c for c in ['title', 'abstract', 'methods_text', 'results_text'] if c in cols]
    output_rows = defaultdict(list)
    years_filter = years_to_write if years_to_write is not None else None

    meta_values = {name: values_for_column(cols, name, n) for name in OUTPUT_META_COLUMNS}

    for i in range(n):
        year = meta_values['year'][i]
        year_key = None if year is None else int(year)
        if years_filter is not None and year_key not in years_filter:
            continue
        if required_text_columns and not row_has_any_value(cols, required_text_columns, i):
            continue

        row = {}
        for name in OUTPUT_META_COLUMNS:
            if name == 'year' and 'year' not in schema.names:
                continue
            row[name] = meta_values[name][i]
        row['has_pmc_full_text'] = int(row_has_any_value(cols, [c for c in FULL_TEXT_MARKER_COLUMNS if c in cols], i))
        row['blinded_text'] = build_blinded_text(cols, i)
        output_rows[year_key].append(row)

    return {
        year: pa.Table.from_pylist(rows, schema=schema)
        for year, rows in output_rows.items()
        if rows
    }


class YearPartitionWriter:
    def __init__(self, out_path, schema, rows_per_file):
        self.out_path = out_path
        self.schema = schema
        self.rows_per_file = max(1, int(rows_per_file))
        self.tmp_tag = f'tmp-{os.getpid()}'
        self.writers = {}
        self.part_rows = defaultdict(int)
        self.part_index = defaultdict(int)
        self.total_rows = defaultdict(int)
        self.tmp_dirs = {}

    def _tmp_dir(self, year):
        if year not in self.tmp_dirs:
            tmp_dir = self.out_path / f'.year={hive_year_value(year)}.{self.tmp_tag}'
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)
            tmp_dir.mkdir(parents=True, exist_ok=True)
            self.tmp_dirs[year] = tmp_dir
        return self.tmp_dirs[year]

    def _part_path(self, year):
        return self._tmp_dir(year) / f'part-{self.part_index[year]:05d}.parquet'

    def _open_writer(self, year):
        writer = self.writers.get(year)
        if writer is None:
            writer = pq.ParquetWriter(self._part_path(year), self.schema, compression='zstd')
            self.writers[year] = writer
        return writer

    def write(self, year, table):
        if table.num_rows == 0:
            return
        writer = self._open_writer(year)
        writer.write_table(table)
        self.part_rows[year] += table.num_rows
        self.total_rows[year] += table.num_rows
        if self.part_rows[year] >= self.rows_per_file:
            writer.close()
            self.writers.pop(year, None)
            self.part_index[year] += 1
            self.part_rows[year] = 0

    def close(self):
        for writer in list(self.writers.values()):
            writer.close()
        self.writers.clear()

        for year, tmp_dir in self.tmp_dirs.items():
            final_dir = partition_dir(self.out_path, year)
            if final_dir.exists():
                shutil.rmtree(final_dir)
            tmp_dir.rename(final_dir)


def run_pyarrow(args, out_path, partitions, where_clause):
    dataset = open_pyarrow_dataset(args.master)
    schema_names = set(dataset.schema.names)
    include_year = not (args.partition_by_year and out_path.suffix != '.parquet')
    schema = output_schema(dataset.schema, include_year=include_year)

    existing_columns = set(dataset.schema.names)
    needed = set(OUTPUT_META_COLUMNS) | {name for name, _ in TEXT_SECTIONS} | set(FULL_TEXT_MARKER_COLUMNS)
    columns = [name for name in dataset.schema.names if name in needed]
    missing = sorted(needed - existing_columns)
    if missing:
        log(f'PyArrow source missing optional columns: {", ".join(missing)}')

    if args.partition_by_year and out_path.suffix != '.parquet':
        years_to_write = partition_years_to_write(out_path, partitions, args.resume)
        if not years_to_write:
            log('All year partitions already complete; no rewrite needed')
            return
        log(f'PyArrow streaming rewrite will write {len(years_to_write):,}/{len(partitions):,} year partitions')
        for year in years_to_write:
            stale = partition_dir(out_path, year)
            if stale.exists() and parquet_row_count(stale, strict=False) <= 0:
                shutil.rmtree(stale)
        filter_expr = pyarrow_filter_expr(schema_names, args.full_text_only, years_to_write)
        scanner = dataset.scanner(columns=columns, filter=filter_expr, batch_size=int(args.batch_size))
        writer = YearPartitionWriter(out_path, schema, args.rows_per_file)
        total = 0
        try:
            for batch in scanner.to_batches():
                for year, table in batch_to_tables_by_year(batch, schema, args.full_text_only, years_to_write).items():
                    writer.write(year, table)
                    total += table.num_rows
                    if total % int(args.progress_every_rows) < table.num_rows:
                        log(f'PyArrow streamed {total:,} rows')
        finally:
            writer.close()
        log(f'PyArrow streamed rows written: {total:,}')
    else:
        filter_expr = pyarrow_filter_expr(schema_names, args.full_text_only, None)
        scanner = dataset.scanner(columns=columns, filter=filter_expr, batch_size=int(args.batch_size))
        target_path = out_path if out_path.suffix == '.parquet' else out_path / 'part-00000.parquet'
        writer = None
        total = 0
        try:
            for batch in scanner.to_batches():
                rows = []
                tables = batch_to_tables_by_year(batch, schema, args.full_text_only, None)
                for table in tables.values():
                    rows.extend(table.to_pylist())
                if not rows:
                    continue
                table = pa.Table.from_pylist(rows, schema=schema)
                if writer is None:
                    writer = pq.ParquetWriter(target_path, schema, compression='zstd')
                writer.write_table(table)
                total += table.num_rows
                if total % int(args.progress_every_rows) < table.num_rows:
                    log(f'PyArrow streamed {total:,} rows')
        finally:
            if writer:
                writer.close()
        log(f'PyArrow streamed rows written: {total:,}')


def run_duckdb(args, out_path, where_clause):
    con = connect_duckdb(args)
    try:
        con.execute(f"CREATE VIEW m AS SELECT * FROM read_parquet('{parquet_source(args.master)}', union_by_name=true)")
        copy_options = 'FORMAT PARQUET, COMPRESSION ZSTD'
        full_text_expr = full_text_expr_sql()
        if args.partition_by_year and out_path.suffix != '.parquet':
            partitions = con.execute(f'''
                SELECT year, count(*) AS n
                FROM m
                WHERE {where_clause}
                GROUP BY year
                ORDER BY year NULLS LAST
            ''').fetchall()
            log(f'DuckDB COPY for blinded text started: {len(partitions):,} year partitions')
            for i, (year, rows) in enumerate(partitions, start=1):
                part_dir = out_path / f'year={hive_year_value(year)}'
                part_dir.mkdir(parents=True, exist_ok=True)
                part_file = part_dir / 'part-00000.parquet'
                log(f'Writing year partition {i:,}/{len(partitions):,}: year={year}; rows={rows:,}')
                con.execute(f'''
                    COPY (
                        {blinded_select(full_text_expr, include_year=False)}
                        WHERE {where_clause} AND {year_filter(year)}
                    ) TO '{sql_string(part_file)}' ({copy_options})
                ''')
        else:
            target_path = out_path
            if out_path.suffix != '.parquet':
                target_path = out_path / 'part-00000.parquet'
            log('DuckDB COPY for blinded text started')
            con.execute(f'''
                COPY (
                    {blinded_select(full_text_expr, include_year=True)}
                    WHERE {where_clause}
                ) TO '{sql_string(target_path)}' ({copy_options})
            ''')
    finally:
        con.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--master', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--temp_dir', default=None)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--memory_limit', default=None)
    ap.add_argument('--partition_by_year', action='store_true')
    ap.add_argument('--full_text_only', action='store_true')
    ap.add_argument('--engine', choices=['pyarrow', 'duckdb'], default='pyarrow')
    ap.add_argument('--batch_size', type=int, default=1000)
    ap.add_argument('--rows_per_file', type=int, default=100000)
    ap.add_argument('--progress_every_rows', type=int, default=250000)
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()

    out_path = Path(args.out)
    directory_output = out_path.suffix != '.parquet'
    if directory_output:
        if out_path.exists() and not args.resume:
            shutil.rmtree(out_path)
        out_path.mkdir(parents=True, exist_ok=True)
        clean_stale_temp_dirs(out_path)
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)

    where_clause = where_clause_sql(args.full_text_only)
    log(f'Reading master table: {args.master}')
    input_rows = parquet_row_count(args.master)
    log(f'master input rows: {input_rows:,}')
    log(f'Writing blinded text: {args.out}')
    memory_label = args.memory_limit or 'default'
    log(
        f'Engine={args.engine}; DuckDB threads={args.threads}; memory_limit={memory_label}; '
        f'partition_by_year={bool(args.partition_by_year)}; full_text_only={bool(args.full_text_only)}; '
        f'resume={bool(args.resume)}'
    )

    if args.engine == 'pyarrow' and args.partition_by_year and directory_output:
        partitions = duckdb_partitions(args, where_clause)
        log(f'DuckDB count for blinded text found {len(partitions):,} year partitions')
        run_pyarrow(args, out_path, partitions, where_clause)
    elif args.engine == 'pyarrow':
        run_pyarrow(args, out_path, [], where_clause)
    else:
        run_duckdb(args, out_path, where_clause)

    log('Blinded text write finished')
    n = parquet_row_count(out_path)
    log(f'blinded_text rows: {n:,}')


if __name__ == '__main__':
    main()
