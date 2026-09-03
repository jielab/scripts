#!/usr/bin/env python3


# 🚩 Imports, paths, and shared inputs
import argparse
import re
import shutil
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from datetime import datetime

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

SAMPLE_RX = re.compile(r'\b(?:n\s*=\s*|N\s*=\s*|sample(?: size)?(?: of)?\s+|participants?\s*[:=]?\s*)([0-9][0-9,]{1,})\b')
MAX_SAMPLE_SIZE = 1_000_000_000

FLAGS = {
    'has_randomized': ['randomized', 'randomised', 'randomly assigned'],
    'has_prospective': ['prospective'],
    'has_cohort': ['cohort'],
    'has_meta_analysis': ['meta-analysis', 'meta analysis', 'systematic review'],
    'has_gwas': ['genome-wide association', 'genome wide association', 'gwas'],
    'has_mr': ['mendelian randomization', 'mendelian randomisation'],
    'has_external_validation': ['external validation', 'validated in an independent', 'independent validation'],
    'has_sensitivity_analysis': ['sensitivity analysis', 'sensitivity analyses'],
    'has_multiple_testing': ['multiple testing', 'false discovery rate', 'bonferroni', 'fdr'],
    'has_data_availability': ['data availability', 'data are available', 'data is available'],
    'has_code_availability': ['code availability', 'github', 'source code', 'code is available'],
    'has_trial_registration': ['trial registration', 'clinicaltrials.gov', 'isrctn'],
    'has_ethics': ['ethics', 'institutional review board', 'informed consent', 'ethical approval'],

    # Paper-score dimensions for weak labels / model audit. These are transparent,
    # imperfect proxies and should be replaced or calibrated by human/LLM labels
    # whenever possible. They intentionally use blinded text only.
    'has_method_innovation': [
        'novel method', 'new method', 'new algorithm', 'new framework', 'new model',
        'we developed', 'we propose', 'we introduce', 'methodological innovation',
        'machine learning', 'deep learning', 'artificial intelligence', 'transformer'
    ],
    'has_new_discovery_language': [
        'we discovered', 'we identified', 'we found that', 'reveals', 'uncovered',
        'previously unknown', 'first evidence', 'first report', 'new insight'
    ],
    'has_large_resource_language': [
        'biobank', 'registry', 'consortium', 'multi-center', 'multicenter',
        'nationwide', 'population-based', 'large-scale', 'atlas', 'database',
        'whole-genome', 'whole genome', 'single-cell', 'multi-omics', 'multiomics'
    ],
    'has_replication': ['replication cohort', 'replicated in', 'replicate our findings', 'validation cohort'],
    'has_open_science_language': ['open access data', 'publicly available', 'github', 'zenodo', 'figshare', 'osf'],
}

FEATURE_SCHEMA = pa.schema([
    ('pmid', pa.string()),
    ('pmcid', pa.string()),
    ('doi', pa.string()),
    ('year', pa.int64()),
    ('word_count', pa.int64()),
    ('sample_size_max', pa.int64()),
    *[(name, pa.int64()) for name in FLAGS],
    ('quality_proxy_0_20', pa.int64()),
])

def log(message):
    now = datetime.now().isoformat(timespec='seconds')
    print(f'[{now}] {message}', flush=True)

def parquet_source(path):
    p = Path(path)
    if p.is_dir() or p.suffix != '.parquet':
        return p
    return p

def open_parquet_dataset(path):
    p = parquet_source(path)
    if Path(p).is_dir():
        return ds.dataset(p, format='parquet', partitioning='hive')
    return ds.dataset(p, format='parquet')

def parquet_row_count(path):
    p = Path(path)
    if p.is_file() and p.suffix == '.parquet':
        return int(pq.ParquetFile(p).metadata.num_rows)
    if not p.exists():
        return 0
    total = 0
    for file_path in p.rglob('*.parquet'):
        total += int(pq.ParquetFile(file_path).metadata.num_rows)
    return total

def safe_sample_size(raw):
    digits = raw.replace(',', '')
    if len(digits) > 18:
        return None
    try:
        value = int(digits)
    except ValueError:
        return None
    if value > MAX_SAMPLE_SIZE:
        return None
    return value

def max_sample_size(text):
    vals = []
    for m in SAMPLE_RX.finditer(text[:20000]):
        value = safe_sample_size(m.group(1))
        if value is not None:
            vals.append(value)
    return max(vals) if vals else None

def quality_proxy(row):
    score = 0
    positive = [
        'has_randomized', 'has_prospective', 'has_cohort', 'has_meta_analysis',
        'has_gwas', 'has_mr', 'has_external_validation', 'has_sensitivity_analysis',
        'has_multiple_testing', 'has_data_availability', 'has_code_availability',
        'has_trial_registration', 'has_ethics',
        'has_method_innovation', 'has_new_discovery_language',
        'has_large_resource_language', 'has_replication', 'has_open_science_language'
    ]
    score += sum(1 for k in positive if row.get(k))
    n = row.get('sample_size_max')
    if n:
        if n >= 1000: score += 2
        elif n >= 100: score += 1
    wc = row.get('word_count') or 0
    if wc >= 1000: score += 2
    elif wc >= 300: score += 1
    return min(score, 20)

def process_batch(batch):
    cols = batch.to_pydict()
    out = []
    n = len(cols['pmid'])
    for i in range(n):
        text = cols.get('blinded_text', [''])[i] or ''
        low = text.lower()
        row = {
            'pmid': cols.get('pmid', [None]*n)[i],
            'pmcid': cols.get('pmcid', [None]*n)[i],
            'doi': cols.get('doi', [None]*n)[i],
            'year': cols.get('year', [None]*n)[i],
            'word_count': len(re.findall(r'\w+', text)),
            'sample_size_max': max_sample_size(text),
        }
        for k, pats in FLAGS.items():
            row[k] = int(any(p in low for p in pats))
        row['quality_proxy_0_20'] = quality_proxy(row)
        out.append(row)
    return pa.Table.from_pylist(out, schema=FEATURE_SCHEMA)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--blinded', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--batch_size', type=int, default=5000)
    ap.add_argument('--rows_per_file', type=int, default=250000)
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--progress_every_rows', type=int, default=250000)
    args = ap.parse_args()

    out_path = Path(args.out)
    directory_output = out_path.suffix != '.parquet'
    if directory_output:
        if out_path.exists():
            shutil.rmtree(out_path)
        out_path.mkdir(parents=True, exist_ok=True)
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)
    input_rows = parquet_row_count(args.blinded)
    log(f'Reading blinded text: {args.blinded}')
    log(f'blinded_text input rows: {input_rows:,}')
    log(
        f'Feature settings: workers={max(1, int(args.workers))}; '
        f'batch_size={args.batch_size:,}; rows_per_file={args.rows_per_file:,}; '
        f'max_sample_size={MAX_SAMPLE_SIZE:,}'
    )
    dataset = open_parquet_dataset(args.blinded)
    writer = None
    total = 0
    part_index = 0
    part_rows = 0
    workers = max(1, int(args.workers))
    pending = set()
    started = time.monotonic()

    def write_features(table):
        nonlocal writer, total, part_index, part_rows
        if writer is None:
            if directory_output:
                part_path = out_path / f'part-{part_index:05d}.parquet'
            else:
                part_path = out_path
            writer = pq.ParquetWriter(part_path, table.schema, compression='zstd')
        writer.write_table(table)
        total += table.num_rows
        part_rows += table.num_rows
        if directory_output and part_rows >= args.rows_per_file:
            writer.close()
            writer = None
            part_index += 1
            part_rows = 0
        progress_every = max(1, int(args.progress_every_rows))
        if total % progress_every < table.num_rows:
            elapsed = max(time.monotonic() - started, 1.0)
            rows_per_min = total / elapsed * 60.0
            remaining = max(input_rows - total, 0)
            eta_min = remaining / rows_per_min if rows_per_min > 0 else 0.0
            log(f'processed {total:,}/{input_rows:,} rows; rate={rows_per_min:,.0f} rows/min; eta={eta_min:.1f} min')

    if workers == 1:
        for batch in dataset.to_batches(batch_size=args.batch_size):
            write_features(process_batch(batch))
    else:
        log(f'Feature extraction workers: {workers}')
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for batch in dataset.to_batches(batch_size=args.batch_size):
                pending.add(executor.submit(process_batch, batch))
                if len(pending) >= workers * 2:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        write_features(future.result())
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    write_features(future.result())
    if writer:
        writer.close()
    log(f'article_features rows: {total:,}')

if __name__ == '__main__':
    main()
