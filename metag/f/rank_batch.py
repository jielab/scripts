"""Reuse the existing resident-index engine without changing ranked evidence."""
from pathlib import Path
import tempfile
from contextlib import nullcontext
from common import atomic_json
from ordinal_query import query_many
from rank_votes import MateVotes


def run_batch(cfg, samples, snapshot, line, replace, run_sample, query_only_dir=None):
    pending = []
    for sample in samples:
        dest = Path(cfg['output_dir']) / sample['dataset'] / sample['sample']
        if dest.exists() and not replace:
            # Use the original complete-result validation, never existence alone.
            run_sample(cfg, sample, snapshot, line, replace)
        else:
            pending.append(sample)
    if not pending:
        if query_only_dir is not None:
            Path(query_only_dir).mkdir(parents=True, exist_ok=False)
            atomic_json(Path(query_only_dir)/'complete.json', dict(release=cfg['rank_release'], samples=samples, pending=[], queries=[]))
        return
    if query_only_dir is None and (len(pending) == 1 or not cfg.get('batch_query', True)):
        for sample in pending:
            run_sample(cfg, sample, snapshot, line, replace)
        return
    root = Path(cfg['work_dir']) / 'tmp' / 'ranked-batches'
    root.mkdir(parents=True, exist_ok=True)
    if query_only_dir is not None:
        Path(query_only_dir).mkdir(parents=True, exist_ok=False)
    context = (nullcontext(str(query_only_dir)) if query_only_dir is not None
               else tempfile.TemporaryDirectory(prefix='query-', dir=root))
    with context as temporary:
        work = Path(temporary)
        review = MateVotes(line)
        queries = query_many(cfg, pending, snapshot, work, taxonomy=line.tax, observe=review.review)
        count = sum(q['fragments'] for q in queries)
        mates = work / 'all-mate-taxids.bin'
        review.finish(mates, count)
        with mates.open('rb') as source:
            for sample, q in zip(pending, queries, strict=True):
                if q['sample_id'] != sample['sample_id']:
                    raise RuntimeError('Resident batch sample order mismatch')
                source.seek((q['first_ordinal'] - 1) * 8)
                target = work / (str(q['first_ordinal']) + '.mate-taxids.bin')
                remaining = q['fragments'] * 8
                with target.open('wb') as output:
                    while remaining:
                        block = source.read(min(8 << 20, remaining))
                        if not block:
                            raise RuntimeError('Incomplete batch mate-vote slice')
                        output.write(block)
                        remaining -= len(block)
                q['mate_taxids'] = str(target)
                if query_only_dir is None:
                    run_sample(cfg, sample, snapshot, line, replace, prepared_query=q)
        if query_only_dir is not None:
            atomic_json(work/'complete.json', dict(release=cfg['rank_release'], samples=samples, pending=pending, queries=queries))
