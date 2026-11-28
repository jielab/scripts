"""Full-library read-pair classification, competitive target evidence and atomic FASTQ publication."""
import argparse
from collections import Counter
from contextlib import contextmanager
import csv
from functools import lru_cache
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

from common import CODE_ROOT, SHARED_RUNTIME, atomic_json, default_cores, fingerprint, load_config, locked, resolve_mode, sha256
from compressed_output import BufferedWriters, gzip_text
from fastq_stream import fragments, new_digest, update_digest
from index_selection import select_index
from native_union import executable_path
from ordinal_query import query_once
from partitioned_build import database_snapshot
from project_cache import batch_plan
from query import discover_samples
from rank_progress import SampleProgress
from rank_votes import MateVotes
from rank_evidence import Lineages, choose_panel, map_all, reconcile
from scheduler import write_cmd, submit_cmd
from storage_layout import new_run
from ranked_storage import publish_directory, recover_directory, remove_entry
from species_split import taxids
from stage_alignment import _mapper
from taxonomy import Taxonomy, safe_name

FORMAT = 'metag-ranked-v2'


def stamp(path):
    path = Path(path).resolve(strict=True)
    st = path.stat()
    return [str(path), st.st_size, st.st_mtime_ns]


def tsv(path, header, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def release(cfg):
    """Record the installed code identity. Workers check it before opening input data."""
    paths = list(CODE_ROOT.glob('*.json')) + list(CODE_ROOT.glob('*.sh'))
    paths += [p for p in (CODE_ROOT/'f').rglob('*') if p.is_file()
              and '__pycache__' not in p.parts and p.suffix != '.pyc']
    paths += [p for p in SHARED_RUNTIME.iterdir() if p.is_file() and p.suffix in ('.py', '.sh', '.yml', '.json')]
    items = [(str(p.relative_to(CODE_ROOT.parent)), sha256(p)) for p in sorted(paths)]
    key = fingerprint(items)
    # A reporting-only change should not force a new multi-terabyte query.
    report_files = {'rank_paper.py', 'paper_layout.py', 'paper_summary.py', 'metag.paper.R'}
    cfg['rank_method'] = fingerprint([(path, digest) for path, digest in items
                          if path.startswith('metag/f/') and Path(path).name not in report_files])
    return key, CODE_ROOT


def effective_config(args):
    mode = resolve_mode(args.submit)
    os.environ['submit'] = os.environ['SUBMIT'] = mode
    cfg = select_index(load_config(args.config, args.root), 'full')
    cfg.update(output_dir=str(Path(cfg['index_base_output_dir'])/'ranked'),
               work_dir=str(Path(cfg['index_base_work_dir'])/'main/ranked'),
               group_library_dir=str(Path(cfg['library_dir'])/'group/optimized'),
               classifier_engine='kraken', union_memory_mapping=cfg.get('rank_union_memory_mapping', False))
    cores = args.cores or default_cores(mode)
    cfg.update(classify_cores=cores, classify_threads=args.threads or max(1, cores-2),
               classify_queue=args.queue or cfg.get('rank_classify_queue', cfg['build_queue']),
               classify_memory_gib=args.mem_gib or cfg.get('rank_classify_memory_gib', 2048),
               classify_hours=args.hours or cfg['classify_hours'])
    opts = cfg.get('stage_alignment') or {}
    cfg.update(rank_min_windows=args.min_windows if args.min_windows is not None else cfg.get('rank_min_windows', 3),
               rank_max_refs=args.max_refs if args.max_refs is not None else opts.get('max_references', 60),
               rank_panel_gib=args.panel_gib if args.panel_gib is not None else opts.get('max_reference_gib', 24))
    raw_taxids = (args.panel_taxids.split(',') if args.panel_taxids else cfg.get('rank_panel_taxids', []))
    cfg['rank_panel_taxids'] = [positive(str(t)) for t in raw_taxids]
    if args.confidence is not None:
        cfg['confidence'] = args.confidence
    table = (Path(args.representatives).expanduser().resolve() if args.representatives else
             Path(cfg['root']) / (opts.get('representatives') or cfg.get('stage_alignment_representatives')
                                  or Path(cfg['group_library_dir']) / 'references.tsv'))
    if not table.is_file():
        raise FileNotFoundError(f'Missing ranked reference catalog: {table}. See --representatives.')
    cfg['stage_alignment'] = dict(opts, representatives=str(table.resolve()))
    if (cfg['classify_threads'] > cores or cfg['rank_min_windows'] < 1 or cfg['rank_max_refs'] < 1
            or not math.isfinite(cfg['rank_panel_gib']) or cfg['rank_panel_gib'] <= 0 or not 0 <= cfg['confidence'] <= 1):
        raise ValueError('Invalid thread limit, confidence, reference budget or minimum windows')
    return cfg


def assignments(cfg, line, query, sparse, eligible, seeds):
    def hits():
        with sparse.open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                yield int(row[0]), (int(row[1]), row[2])
    stream = hits()
    next_hit = next(stream, None)
    decide = lru_cache(maxsize=131072)(lambda k, s: reconcile(
        k, s, line, eligible, seeds, cfg['rank_min_windows']))
    count = query['fragments']
    calls = zip(taxids(query['taxids'], count), taxids(query['mate_taxids'], count), strict=True)
    for i, (k, mate) in enumerate(calls, 1):
        hit = None
        if next_hit and i == next_hit[0]:
            hit = next_hit[1]
            next_hit = next(stream, None)
        t, reason = decide(k, hit)
        if not hit and line.tax.rank.get(k) == 'species':
            if mate != k:
                t = line.common([t, mate]) if mate else 0
                reason = 'mate_vote_downgrade'
            elif reason == 'classifier_only':
                reason = 'both_mates_kraken'
        yield i, k, mate, hit, t, reason
    if next_hit:
        raise RuntimeError('Unconsumed or unordered panel evidence')


def emit(cfg, sample, line, query, sparse, eligible, seeds, out, progress):
    @lru_cache(maxsize=None)
    def label(t):
        rank = line.tax.rank.get(t, 'unclassified')
        group = line.group(t)
        return group, rank, line.tax.name.get(t, 'unclassified'), f'bins/{group}/{safe_name(rank)}/{t}'
    bins = Counter()
    reasons = Counter()
    taxon_evidence = Counter()
    transfers = Counter()
    hashes = {}
    digest = new_digest()
    total = 0
    writers = BufferedWriters(out, threads=min(
        cfg['classify_threads'], 8), compresslevel=cfg['split_gzip_level'])
    try:
        with gzip_text(out/'assignments.tsv.gz', threads=1) as audit:
            audit.write('ordinal\tkraken_taxid\tmate_review_taxid\tpanel_taxid\tpanel_evidence\tfinal_taxid\tevidence\n')
            raw = fragments(sample['r1'], sample['r2'])
            for (a, b), (i, k, mate, hit, t, reason) in zip(raw, assignments(cfg, line, query, sparse, eligible, seeds), strict=True):
                group, rank, name, base = label(t)
                update_digest(digest, a, b)
                if t not in hashes:
                    hashes[t] = new_digest()
                update_digest(hashes[t], a, b)
                writers.write(base+'/R1.fastq.gz', a)
                writers.write(base+'/R2.fastq.gz', b)
                audit.write(
                    f'{i}\t{k}\t{mate}\t{hit[0] if hit else 0}\t{hit[1] if hit else "none"}\t{t}\t{reason}\n')
                bins[t] += 1
                reasons[reason] += 1
                taxon_evidence[t, reason] += 1
                transfers[(label(k)[0], group)] += 1
                total = i
                if i % 100000 == 0:
                    progress.counted(i)
        if total != query['fragments'] or digest.hexdigest() != query['input_digest']:
            raise RuntimeError('FASTQ changed between classification and split')
    except BaseException:
        writers.close(flush=False)
        raise
    progress.counted(total)
    progress.stage('flush_bins', bins=len(bins))
    writers.close()
    counts = []
    for t, n in bins.most_common():
        group, rank, name, base = label(t)
        sizes = [(out/base/f'R{i}.fastq.gz').stat().st_size for i in (1, 2)]
        counts.append([t, name, group, rank, n, base, *sizes, hashes[t].hexdigest()])
    tsv(out/'taxa.tsv', ['taxid', 'name', 'group', 'rank', 'pairs',
        'bin', 'r1_bytes', 'r2_bytes', 'pairs_sha256'], counts)
    tsv(out/'taxon_evidence.tsv', ['final_taxid', 'evidence', 'pairs'], [(*k, n)
        for k, n in sorted(taxon_evidence.items())])
    tsv(out/'evidence.tsv', ['evidence', 'pairs'], sorted(reasons.items()))
    tsv(out/'transitions.tsv', ['kraken_group', 'final_group', 'pairs'], [(*k, n)
        for k, n in sorted(transfers.items())])
    return dict(pairs=total, exact_species_pairs=sum(r[4] for r in counts if r[3] == 'species'),
                unclassified_pairs=bins[0], input_digest=digest.hexdigest())


@contextmanager
def work_directory(gen, cfg, sample):
    with tempfile.TemporaryDirectory(dir=gen, prefix='work-') as tmp:
        work = Path(tmp)
        try:
            yield work
        except BaseException:
            logs = Path(cfg['work_dir'])/'logs/failed'/sample['dataset']/sample['sample']/gen.name
            logs.mkdir(parents=True, exist_ok=True)
            for path in work.glob('*.log'):
                shutil.copyfile(path, logs/path.name)
            print('Failure logs: '+str(logs), file=sys.stderr, flush=True)
            raise


def run_sample(cfg, sample, snapshot, line, replace, prepared_query=None):
    dest = Path(cfg['output_dir'])/sample['dataset']/sample['sample']
    gens = dest.parent/('.'+dest.name+'.runs')
    opts = cfg.get('stage_alignment') or {}
    native = executable_path(cfg)
    mapper = _mapper(cfg, opts)
    context = dict(code=cfg['rank_method'], database=snapshot['identity'], raw=[stamp(sample[k]) for k in ('r1', 'r2')],
                   programs=[stamp(native), stamp(mapper)], catalog=[stamp(opts['representatives'])],
                   settings={k: cfg[k] for k in ('rank_panel_gib', 'rank_max_refs', 'rank_min_windows', 'rank_panel_taxids', 'confidence',
                                                 'minimum_hit_groups', 'minimum_base_quality')}, reference_raw=cfg['reference_raw'],
                   metadata_dir=cfg['metadata_dir'], metadata=[stamp(p) for p in sorted((Path(cfg['metadata_dir'])/'summaries').glob('*')) if p.is_file()])
    with locked(gens/'write.lock'):
        with locked(gens/'publish.lock'):
            recover_directory(dest)
        if dest.is_symlink():
            raise RuntimeError('Expected a regular ranked sample directory: '+str(dest))
        if dest.exists():
            receipt = dest / 'summary.json'
            record = json.loads(receipt.read_text()) if receipt.is_file() else {}
            if record.get('format') != FORMAT or record.get('sample') != sample['sample_id']:
                raise RuntimeError('Directory already in use: ' + str(dest))
        if dest.exists() and not replace:
            previous = dict(record.get('request') or {})
            if (previous != context
                    or any(stamp(s[0]) != s for s in record['reference_states'])):
                raise RuntimeError(
                    'Inputs, references or settings changed; use --replace TRUE')
            for name, size in record['file_sizes'].items():
                if (dest/name).stat().st_size != size:
                    raise RuntimeError('Incomplete output: '+name+'; use --replace TRUE')
            print('Reuse '+str(dest), flush=True)
            return
        gen = Path(tempfile.mkdtemp(prefix='run-', dir=gens))
        try:
            with work_directory(gen, cfg, sample) as work, SampleProgress(sample['sample_id']) as progress:
                progress.stage('query_input', prepared=prepared_query is not None)
                if prepared_query is None:
                    review = MateVotes(line)
                    q = query_once(cfg, sample, snapshot, work, taxonomy=line.tax, observe=review.review)
                    q['mate_taxids'] = str(work/'mate_taxids.bin')
                    review.finish(q['mate_taxids'], q['fragments'])
                else:
                    from ordinal_query import input_states
                    q = dict(prepared_query)
                    if q['sample_id'] != sample['sample_id'] or q['input_states'] != input_states(sample):
                        raise RuntimeError('Prepared query does not match unchanged sample inputs')
                    shutil.copyfile(q['log_path'], work/'query.log')
                count = q['fragments']
                progress.stage('select_panel', expected=count)
                counts = Counter(taxids(q['taxids'], count))
                refs, eligible, omitted, table = choose_panel(cfg, line, counts, work)
                states = [stamp(r['resolved_path']) for r in refs]
                sparse, seeds, mapping = map_all(cfg, sample, refs, zip(taxids(q['taxids'], count), taxids(
                    q['mate_taxids'], count), strict=True), line, work, count, q['input_digest'], progress)
                progress.stage('write_bins', expected=count)
                stats = emit(cfg, sample, line, q, sparse, eligible, seeds, gen, progress)
                progress.stage('finish_tables')
                if states != [stamp(r['resolved_path']) for r in refs]:
                    raise RuntimeError('Selected FASTA changed during run')
                tsv(gen/'panel.tsv', ['taxid', 'name', 'role', 'bases', 'kraken_pairs', 'near_relative_and_background', 'agreeing_nuclear_windows_capped', 'reference'],
                    [[r['species_taxid'], r['name'], r['role'], r['bases'], counts[r['species_taxid']], int(r['species_taxid'] in eligible),
                      seeds.get(r['species_taxid'], 0), r['resolved_path']] for r in refs])
                tsv(gen/'panel_omitted.tsv', ['taxid', 'kraken_pairs', 'reason'],
                    [[t, counts[t], reason] for t, reason in sorted(omitted.items())])
                for name in ('query.log', 'mapping.log'):
                    shutil.copyfile(work/name, gen/name)
                sizes = {str(p.relative_to(gen)): p.stat().st_size for p in gen.iterdir() if p.is_file()}
            with (gen/'taxa.tsv').open() as h:
                for row in csv.DictReader(h, delimiter='\t'):
                    for m in (1, 2):
                        sizes[row['bin']+f'/R{m}.fastq.gz'] = int(row[f'r{m}_bytes'])
            atomic_json(gen/'summary.json', dict(format=FORMAT, release=cfg['rank_release'], sample=sample['sample_id'], request=context, taxonomy=snapshot['taxonomy_dir'],
                                                 reference_states=states, panel_catalog=table,
                                                 phase_seconds=progress.timing, mapping=mapping, kraken_command=q['command'], file_sizes=sizes, **stats))
            with locked(gens/'publish.lock'):
                publish_directory(gen, dest)
            print('Published '+str(dest), flush=True)
        finally:
            remove_entry(gen)


def positive(s):
    n = int(s)
    if n < 1:
        raise argparse.ArgumentTypeError('Expected a positive integer')
    return n


def main():
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        job = json.loads(Path(sys.argv[2]).read_text())
        cfg = job['config']
        if release(dict(cfg))[0] != cfg['rank_release']:
            raise RuntimeError('Installed scripts changed after planning. Submit this job again.')
        snap = database_snapshot(cfg)
        line = Lineages(Taxonomy(snap['taxonomy_dir']))
        from rank_batch import run_batch
        stage = job['stage']
        if stage == 'query':
            run_batch(cfg, job['samples'], snap, line, job['replace'], run_sample,
                      query_only_dir=job['prepared_dir'])
        elif stage == 'finish':
            prepared = Path(job['prepared_dir'])
            receipt = json.loads((prepared/'complete.json').read_text())
            if receipt['release'] != cfg['rank_release'] or receipt['samples'] != job['samples']:
                raise RuntimeError('Prepared batch does not match this job')
            for sample, query in zip(receipt['pending'], receipt['queries'], strict=True):
                run_sample(cfg, sample, snap, line, job['replace'], prepared_query=query)
            # Only remove this job's compact intermediate files after all results publish.
            shutil.rmtree(prepared)
        else:
            run_batch(cfg, job['samples'], snap, line, job['replace'], run_sample)
        return
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'split'))
    p.add_argument('--sample', action='append')
    for name in ('dat-list', 'batch', 'config', 'root', 'representatives', 'panel-taxids', 'queue'):
        p.add_argument('--'+name)
    for name in ('cores', 'threads', 'mem-gib', 'hours', 'batch-size', 'min-windows', 'max-refs'):
        p.add_argument('--'+name, type=positive)
    p.add_argument('--submit', choices=['local', 'hpc', 'auto'])
    p.add_argument('--confidence', type=float)
    p.add_argument('--panel-gib', type=float)
    p.add_argument('--replace', choices=['TRUE', 'FALSE'], default='FALSE')
    p.add_argument('--plan-only', action='store_true')
    a = p.parse_args()
    if a.sample and (a.batch or a.batch_size):
        p.error('--sample excludes --batch and --batch-size')
    cfg = effective_config(a)
    samples = discover_samples(cfg, datasets=a.dat_list, only=a.sample, probe=False)
    if not samples:
        p.error('No complete paired-end samples selected')
    # Discover every complete pair by default. Explicit batch numbers use the requested batch size.
    automatic = not (a.sample or a.batch or a.batch_size)
    if automatic:
        from rank_plan import automatic_batches
        plans = automatic_batches(samples, cfg)
    else:
        plans = batch_plan(cfg, samples, 1 if a.sample else a.batch_size or cfg['sample_batch_size'], a.batch or 'all')
    staged = resolve_mode(a.submit) == 'hpc' and automatic
    print(f'Selected {len(samples)} samples; {len(plans)} resident-index batches; staged={staged}', flush=True)
    key, source = release(cfg)
    cfg.update(rank_release=key, code_root=str(source))
    jobs = Path(cfg['work_dir'])/'jobs'
    jobs.mkdir(parents=True, exist_ok=True)
    folder = new_run(jobs / 'split')
    commands = []
    for dataset, index, members in plans:
        unit = safe_name(members[0]['sample']) if len(members) == 1 else f'batch{index:03d}'
        label = safe_name(dataset)+'.'+unit
        stages = ('query', 'finish') if staged else ('combined',)
        pair = []
        for stage in stages:
            stage_cfg = dict(cfg)
            if stage == 'finish':
                stage_cfg.update(classify_queue=cfg.get('rank_finish_queue', '33'),
                                 classify_memory_gib=cfg.get('rank_finish_memory_gib', 80))
            name = label + ('.' + stage if staged else '')
            payload = dict(config=stage_cfg, samples=members, replace=a.replace == 'TRUE',
                           stage=stage, prepared_dir=str(folder/(label+'.prepared')))
            spec = folder/(name+'.job.json')
            atomic_json(spec, payload)
            command = [sys.executable, '-B', str(source/'f/rank_workflow.py'), '--worker', str(spec)]
            cmd = write_cmd(stage_cfg, folder/(name+'.cmd'), 'classify', command, 'metag.'+name)
            pair.append(cmd)
            print(str(cmd), flush=True)
        commands.append(pair)
    if not a.plan_only:
        for pair in commands:
            row = submit_cmd(cfg, pair[0])
            if staged:
                submit_cmd(cfg, pair[1], dependency=row['job_id'])




if __name__ == '__main__':
    main()
