"""Blind review export with a separate private identity key and real rank provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import uuid

import duckdb
import pandas as pd
from make_pfi_training_dataset import parquet_source, normalize_rank_file, norm_id


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('master', 'blinded', 'features', 'journal_rank', 'out_csv', 'out_xlsx', 'key_out', 'summary_out'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--n', type=int, default=1000)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--min_text_chars', type=int, default=300)
    p.add_argument('--max_scan_rows', type=int, default=200000)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--high_anchor_csv')
    p.add_argument('--low_anchor_csv')
    a = p.parse_args()
    if min(a.n, a.max_scan_rows, a.threads) < 1:
        p.error('n, max_scan_rows and threads must be positive')
    paths = [Path(getattr(a, n)).resolve() for n in ('out_csv','out_xlsx','key_out','summary_out')]
    if len(set(paths)) != len(paths):
        raise ValueError('Review outputs and private key must have different paths')
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    rank_out = Path(a.key_out).with_suffix('.rank.parquet')
    rank = normalize_rank_file(a.journal_rank, rank_out)
    con = duckdb.connect()
    con.execute(f'PRAGMA threads={a.threads}')
    con.register('rank_table', rank)
    # Read text only from the pre-blinded corpus. Never substitute master text.
    q = f"""SELECT coalesce(CAST(m.pmid AS VARCHAR),CAST(m.pmcid AS VARCHAR),CAST(m.doi AS VARCHAR)) article_key,
        CAST(m.pmid AS VARCHAR) pmid, CAST(m.pmcid AS VARCHAR) pmcid, CAST(m.doi AS VARCHAR) doi,
        m.journal, b.blinded_text
        FROM read_parquet('{parquet_source(a.master)}',union_by_name=true) m
        JOIN read_parquet('{parquet_source(a.blinded)}',union_by_name=true) b
        ON coalesce(CAST(m.pmid AS VARCHAR),CAST(m.pmcid AS VARCHAR),CAST(m.doi AS VARCHAR)) =
           coalesce(CAST(b.pmid AS VARCHAR),CAST(b.pmcid AS VARCHAR),CAST(b.doi AS VARCHAR))
        WHERE length(b.blinded_text)>={a.min_text_chars}
        LIMIT {a.max_scan_rows}"""
    df = con.execute(q).fetchdf()
    con.close()
    from make_pfi_training_dataset import norm_journal
    df['journal_key'] = df.journal.map(norm_journal)
    df = df.merge(rank[['journal_key']], on='journal_key', how='inner', validate='many_to_one')
    if df.article_key.isna().any() or df.article_key.duplicated().any():
        raise ValueError('Master/blinded join must have one row per article_key')
    anchors = {}
    for path, label in ((a.high_anchor_csv,'high'),(a.low_anchor_csv,'low')):
        if not path:
            continue
        rows = pd.read_csv(path, dtype=str)
        if 'article_key' not in rows:
            from make_pfi_training_dataset import article_key_from
            for col in ('pmid','pmcid','doi'):
                if col not in rows: rows[col] = None
            rows['article_key'] = article_key_from(rows)
        for key in rows.article_key.map(norm_id).dropna():
            if key in anchors and anchors[key] != label:
                raise ValueError('Conflicting high/low anchor: ' + key)
            anchors[key] = label
    missing = set(anchors) - set(df.article_key)
    if missing:
        raise ValueError(f'{len(missing)} anchors absent from bounded eligible blinded corpus; increase --max_scan_rows or correct inputs')
    forced = df[df.article_key.isin(anchors)]
    if len(forced) > a.n:
        raise ValueError('More anchors than requested review rows')
    rest = df[~df.article_key.isin(anchors)]
    chosen = pd.concat([forced, rest.sample(min(a.n-len(forced),len(rest)),random_state=a.seed)])
    if chosen.empty:
        raise ValueError('No blinded texts matched the real rank file')
    chosen = chosen.sample(frac=1,random_state=a.seed).reset_index(drop=True)
    chosen['review_id'] = [uuid.uuid4().hex for _ in range(len(chosen))]
    def redact(row):
        text = str(row.blinded_text)
        for field in ('journal','pmid','pmcid','doi'):
            value = norm_id(row[field])
            if value:
                text = re.sub(r'(?<!\w)' + re.escape(value) + r'(?!\w)', '[redacted]', text, flags=re.I)
        return text
    review = pd.DataFrame({'review_id': chosen.review_id, 'blinded_text': chosen.apply(redact,axis=1),
                           'paper_score': '', 'review_notes': ''})
    key = chosen[['review_id','article_key','pmid','pmcid','doi','journal']].copy()
    key['anchor_type'] = chosen.article_key.map(anchors).fillna('')
    review.to_csv(a.out_csv,index=False)
    review.to_excel(a.out_xlsx,index=False)
    key.to_csv(a.key_out,index=False)
    Path(a.key_out).chmod(0o600)
    summary = dict(schema='pfi-blind-review-v2', n=len(review), bounded_scan_rows=a.max_scan_rows,
                   seed=a.seed, rank_input=str(Path(a.journal_rank).resolve()),
                   rank_sha256=hashlib.sha256(Path(a.journal_rank).read_bytes()).hexdigest(),
                   text_source=str(Path(a.blinded).resolve()), key_private=True,
                   scores='blank_for_review', anchors_private=len(forced))
    Path(a.summary_out).write_text(json.dumps(summary,indent=2)+'\n')
    print(f'Wrote {len(review)} blinded review records and a separate private key')


if __name__ == '__main__':
    main()
