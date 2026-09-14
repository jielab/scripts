#!/usr/bin/env python3


# 🚩 Imports, paths, and shared inputs
import argparse
import io
import multiprocessing as mp
import os
import re
import tarfile
from pathlib import Path

import polars as pl
from lxml import etree
from tqdm import tqdm

def clean_text(x):
    if x is None:
        return None
    if hasattr(x, 'itertext'):
        x = ' '.join(x.itertext())
    x = re.sub(r'\s+', ' ', str(x)).strip()
    return x or None

def xpath_first(root, expr):
    vals = root.xpath(expr)
    if not vals:
        return None
    return clean_text(vals[0])

def xpath_all_text(root, expr):
    vals = []
    for x in root.xpath(expr):
        t = clean_text(x)
        if t:
            vals.append(t)
    return ' '.join(vals) if vals else None

def get_article_id(root, id_type):
    expr = f"//*[local-name()='article-id' and translate(@pub-id-type, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')='{id_type.lower()}']"
    return xpath_first(root, expr)

def get_year(root):
    exprs = [
        "(//*[local-name()='article-meta']//*[local-name()='pub-date']/*[local-name()='year'])[1]",
        "(//*[local-name()='journal-meta']//*[local-name()='pub-date']/*[local-name()='year'])[1]",
    ]
    for e in exprs:
        y = xpath_first(root, e)
        if y and y.isdigit():
            return int(y)
    txt = xpath_all_text(root, "//*[local-name()='pub-date']")
    if txt:
        m = re.search(r'(19|20)\d{2}', txt)
        if m:
            return int(m.group(0))
    return None

def get_license(root):
    vals = root.xpath("//*[local-name()='license']")
    if not vals:
        return None
    lic = vals[0]
    href = lic.get('{http://www.w3.org/1999/xlink}href') or lic.get('href')
    txt = clean_text(lic)
    return href or txt

def get_section_text(root, keywords):
    hits = []
    for sec in root.xpath("//*[local-name()='sec']"):
        title = xpath_first(sec, "./*[local-name()='title']") or ''
        title_l = title.lower()
        if any(k in title_l for k in keywords):
            txt = clean_text(sec)
            if txt:
                hits.append(txt)
    return ' '.join(hits) if hits else None

def parse_xml_bytes(b, max_chars):
    try:
        root = etree.parse(io.BytesIO(b), etree.XMLParser(recover=True, huge_tree=True)).getroot()
    except Exception:
        return None

    year = get_year(root)
    title = xpath_first(root, "(//*[local-name()='article-title'])[1]")
    abstract = xpath_all_text(root, "//*[local-name()='abstract']")
    full_text = xpath_all_text(root, "//*[local-name()='body']")
    if full_text and max_chars > 0:
        full_text = full_text[:max_chars]

    pmcid = get_article_id(root, 'pmc')
    if pmcid and not pmcid.upper().startswith('PMC'):
        pmcid = 'PMC' + pmcid

    doi = get_article_id(root, 'doi')
    if doi:
        doi = doi.lower()

    journal = xpath_first(root, "(//*[local-name()='journal-title'])[1]")
    if not journal:
        journal = xpath_first(root, "(//*[local-name()='journal-id'])[1]")

    return {
        'pmcid': pmcid,
        'pmid': get_article_id(root, 'pmid'),
        'doi': doi,
        'year': year,
        'title': title,
        'abstract': abstract,
        'journal': journal,
        'article_type': root.get('article-type'),
        'license': get_license(root),
        'ref_count': len(root.xpath("//*[local-name()='ref-list']//*[local-name()='ref']")),
        'table_count': len(root.xpath("//*[local-name()='table-wrap']")),
        'fig_count': len(root.xpath("//*[local-name()='fig']")),
        'intro_text': get_section_text(root, ['intro', 'background']),
        'methods_text': get_section_text(root, ['method', 'material', 'statistical analysis']),
        'results_text': get_section_text(root, ['result', 'finding']),
        'discussion_text': get_section_text(root, ['discussion', 'conclusion']),
        'data_availability_text': get_section_text(root, ['data availability', 'availability of data', 'data sharing']),
        'ethics_text': get_section_text(root, ['ethic', 'institutional review', 'consent']),
        'trial_registration_text': get_section_text(root, ['trial registration', 'registration']),
        'full_text': full_text,
    }

def article_schema():
    return {
        'pmcid': pl.Utf8, 'pmid': pl.Utf8, 'doi': pl.Utf8, 'year': pl.Int64,
        'title': pl.Utf8, 'abstract': pl.Utf8, 'journal': pl.Utf8, 'article_type': pl.Utf8,
        'license': pl.Utf8, 'ref_count': pl.Int64, 'table_count': pl.Int64, 'fig_count': pl.Int64,
        'intro_text': pl.Utf8, 'methods_text': pl.Utf8, 'results_text': pl.Utf8,
        'discussion_text': pl.Utf8, 'data_availability_text': pl.Utf8,
        'ethics_text': pl.Utf8, 'trial_registration_text': pl.Utf8, 'full_text': pl.Utf8,
    }

def write_batch(rows, out_file, schema):
    if rows:
        pl.DataFrame(rows, schema_overrides=schema).write_parquet(out_file, compression='zstd')
    else:
        pl.DataFrame({k: [] for k in schema}, schema=schema).write_parquet(out_file, compression='zstd')

def parse_tar(tar_path, out_dir, min_year, max_chars, batch_size):
    tar_path = Path(tar_path)
    out_dir = Path(out_dir)
    done_file = out_dir / f'{tar_path.name}.done'
    part_glob = f'{tar_path.name}.part*.parquet'
    if done_file.exists():
        return str(tar_path), 'skip', 0, 0

    for stale_part in out_dir.glob(part_glob):
        stale_part.unlink()

    rows = []
    part_idx = 0
    row_count = 0
    schema = article_schema()

    def flush_rows(force_empty=False):
        nonlocal rows, part_idx
        if not rows and not force_empty:
            return
        tmp_file = out_dir / f'{tar_path.name}.part{part_idx:05d}.parquet.tmp'
        out_file = out_dir / f'{tar_path.name}.part{part_idx:05d}.parquet'
        write_batch(rows, tmp_file, schema)
        os.replace(tmp_file, out_file)
        rows = []
        part_idx += 1

    try:
        with tarfile.open(tar_path, 'r:gz') as tar:
            for member in tar:
                if not member.isfile() or not member.name.lower().endswith(('.xml', '.nxml')):
                    continue
                f = tar.extractfile(member)
                if f is None:
                    continue
                row = parse_xml_bytes(f.read(), max_chars)
                if row is None:
                    continue
                if row.get('year') is None or row.get('year') < min_year:
                    continue
                rows.append(row)
                row_count += 1
                if len(rows) >= batch_size:
                    flush_rows()
    except Exception as e:
        print(f'ERROR parsing {tar_path}: {e}')
        for stale_part in out_dir.glob(part_glob):
            stale_part.unlink(missing_ok=True)
        return str(tar_path), 'error', row_count, part_idx

    flush_rows(force_empty=(row_count == 0))
    done_file.write_text(f'rows={row_count}\nparts={part_idx}\n', encoding='utf-8')
    return str(tar_path), 'ok', row_count, part_idx

def _worker(args):
    return parse_tar(*args)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in_dir', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--min_year', type=int, default=2000)
    ap.add_argument('--max_chars', type=int, default=50000)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--max_files', type=int, default=0)
    ap.add_argument('--batch_size', type=int, default=500)
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    files = sorted(Path(args.in_dir).glob('*.tar.gz'))
    if args.max_files and args.max_files > 0:
        files = files[:args.max_files]
    print(f'Parsing {len(files)} PMC OA tar.gz files from {args.in_dir}')

    tasks = [(str(f), args.out_dir, args.min_year, args.max_chars, args.batch_size) for f in files]
    counts = {'ok': 0, 'skip': 0, 'error': 0}
    if args.workers <= 1:
        for t in tqdm(tasks):
            _, status, _, _ = parse_tar(*t)
            counts[status] = counts.get(status, 0) + 1
    else:
        with mp.Pool(args.workers, maxtasksperchild=1) as pool:
            for _, status, _, _ in tqdm(pool.imap_unordered(_worker, tasks, chunksize=1), total=len(tasks)):
                counts[status] = counts.get(status, 0) + 1
    print(f"Done: ok={counts.get('ok', 0)} skip={counts.get('skip', 0)} error={counts.get('error', 0)}")

if __name__ == '__main__':
    main()
