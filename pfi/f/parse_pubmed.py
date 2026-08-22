#!/usr/bin/env python3


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import argparse
import gzip
import multiprocessing as mp
import re
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

def find_text(root, path):
    return clean_text(root.find(path))

def get_year(article):
    for path in [
        './/Article/ArticleDate/Year',
        './/JournalIssue/PubDate/Year',
        './/DateCompleted/Year',
        './/DateRevised/Year',
    ]:
        y = find_text(article, path)
        if y and y.isdigit():
            return int(y)
    md = find_text(article, './/JournalIssue/PubDate/MedlineDate')
    if md:
        m = re.search(r'(19|20)\d{2}', md)
        if m:
            return int(m.group(0))
    return None

def get_article_id(article, id_type):
    for x in article.findall('.//ArticleId'):
        if (x.get('IdType') or '').lower() == id_type.lower():
            return clean_text(x)
    return None

def get_doi(article):
    doi = get_article_id(article, 'doi')
    if doi:
        return doi.lower()
    for x in article.findall('.//ELocationID'):
        if (x.get('EIdType') or '').lower() == 'doi':
            v = clean_text(x)
            return v.lower() if v else None
    return None

def get_abstract(article):
    parts = []
    for x in article.findall('.//Abstract/AbstractText'):
        label = x.get('Label')
        txt = clean_text(x)
        if txt:
            parts.append(f'{label}: {txt}' if label else txt)
    return ' '.join(parts) if parts else None

def join_texts(article, path):
    vals = []
    for x in article.findall(path):
        t = clean_text(x)
        if t:
            vals.append(t)
    return '|'.join(vals) if vals else None

def get_authors(article):
    authors = []
    affils = []
    for au in article.findall('.//AuthorList/Author'):
        collective = find_text(au, 'CollectiveName')
        last = find_text(au, 'LastName')
        fore = find_text(au, 'ForeName')
        initials = find_text(au, 'Initials')
        if collective:
            name = collective
        elif last and fore:
            name = f'{fore} {last}'
        elif last and initials:
            name = f'{initials} {last}'
        elif last:
            name = last
        else:
            continue
        authors.append(name)
        for af in au.findall('.//AffiliationInfo/Affiliation'):
            a = clean_text(af)
            if a:
                affils.append(a)
    return {
        'authors': '|'.join(authors) if authors else None,
        'first_author': authors[0] if authors else None,
        'last_author': authors[-1] if authors else None,
        'author_count': len(authors),
        'affiliations': ' || '.join(affils) if affils else None,
    }

def parse_one(xml_gz, out_dir, min_year):
    xml_gz = Path(xml_gz)
    out_dir = Path(out_dir)
    out_file = out_dir / f'{xml_gz.name}.parquet'
    if out_file.exists() and out_file.stat().st_size > 0:
        return str(out_file), 'skip'

    rows = []
    with gzip.open(xml_gz, 'rb') as f:
        context = etree.iterparse(f, events=('end',), tag='PubmedArticle', recover=True, huge_tree=True)
        for _, article in context:
            year = get_year(article)
            if year is not None and year >= min_year:
                pmid = find_text(article, './/MedlineCitation/PMID')
                rows.append({
                    'pmid': pmid,
                    'pmcid': get_article_id(article, 'pmc'),
                    'doi': get_doi(article),
                    'year': year,
                    'title': find_text(article, './/ArticleTitle'),
                    'abstract': get_abstract(article),
                    'journal': find_text(article, './/Journal/Title'),
                    'iso_journal': find_text(article, './/Journal/ISOAbbreviation'),
                    'issn': find_text(article, './/Journal/ISSN'),
                    'publication_types': join_texts(article, './/PublicationTypeList/PublicationType'),
                    'mesh_terms': join_texts(article, './/MeshHeading/DescriptorName'),
                    'keywords': join_texts(article, './/KeywordList/Keyword'),
                    **get_authors(article),
                })
            article.clear()
            parent = article.getparent()
            while parent is not None and article.getprevious() is not None:
                del parent[0]
    schema = {
        'pmid': pl.Utf8, 'pmcid': pl.Utf8, 'doi': pl.Utf8, 'year': pl.Int64,
        'title': pl.Utf8, 'abstract': pl.Utf8, 'journal': pl.Utf8,
        'iso_journal': pl.Utf8, 'issn': pl.Utf8, 'publication_types': pl.Utf8,
        'mesh_terms': pl.Utf8, 'keywords': pl.Utf8, 'authors': pl.Utf8,
        'first_author': pl.Utf8, 'last_author': pl.Utf8, 'author_count': pl.Int64,
        'affiliations': pl.Utf8,
    }
    if rows:
        pl.DataFrame(rows, schema_overrides=schema).write_parquet(out_file, compression='zstd')
    else:
        pl.DataFrame({k: [] for k in schema}, schema=schema).write_parquet(out_file, compression='zstd')
    return str(out_file), 'ok'

def _worker(args):
    return parse_one(*args)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in_dir', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--min_year', type=int, default=2000)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--max_files', type=int, default=0)
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    files = sorted(Path(args.in_dir).glob('pubmed*.xml.gz'))
    if args.max_files and args.max_files > 0:
        files = files[:args.max_files]
    print(f'Parsing {len(files)} PubMed XML files from {args.in_dir}')

    tasks = [(str(f), args.out_dir, args.min_year) for f in files]
    if args.workers <= 1:
        for t in tqdm(tasks):
            print(parse_one(*t))
    else:
        with mp.Pool(args.workers) as pool:
            for out, status in tqdm(pool.imap_unordered(_worker, tasks), total=len(tasks)):
                pass
    print('Done')

if __name__ == '__main__':
    main()
