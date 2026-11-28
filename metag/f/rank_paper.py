"""Report ranked FASTQ outputs; truth and evidence agreement have separate denominators."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import fcntl
import json
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import traceback
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from datetime import datetime, timezone

from common import CODE_ROOT, SHARED_RUNTIME, atomic_json, fingerprint, load_config, locked, resolve_mode, sha256
from fastq_stream import fragments, new_digest, update_digest
from rank_workflow import FORMAT, positive, release
from query import _name
from scheduler import write_cmd, submit_cmd
from storage_layout import new_run
from paper_layout import flat_ranked, paper_root
from paper_summary import render
from ranked_storage import publish_directory, recover_directory, staging
from taxonomy import Taxonomy
from rank_evidence import Lineages

sys.path.insert(0, str(CODE_ROOT / "refine"))
import metag_refine_tables as refine_tables
import seal_audit
sys.path.insert(0,str(CODE_ROOT/"paper"))
import explorer


def rows(path):
    with (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='') as h:
        yield from csv.DictReader(h, delimiter='\t')


def truth_rows(path, sample, single):
    last = 0
    for row in rows(path):
        if 'sample' in row:
            if row['sample'] != sample:
                continue
        elif not single:
            raise ValueError('Multi-sample truth requires a sample column (DATASET/SAMPLE)')
        i, t = int(row['ordinal']), int(row['true_taxid'])
        if i <= last or t < 0:
            raise ValueError(f'Truth must have unique, increasing positive ordinals per sample: {sample} {i}')
        last = i
        yield i, t


class Accuracy:
    def __init__(self, line):
        self.line = line
        self.matrix = Counter()
        self.n = 0

    def add(self, truth, pred):
        line = self.line
        truth = line.tax.resolve(truth)
        self.n += 1
        for rank in ('species', 'genus', 'family', 'group'):
            if rank == 'group':
                a, b = line.group(truth), line.group(pred)
                if a in ('shared', 'unclassified'):
                    a = None
                if b in ('shared', 'unclassified'):
                    b = None
            else:
                a = line.ancestor(truth, {rank}) or None
                b = line.ancestor(pred, {rank}) or None
            # A nonzero truth known only to a higher rank cannot adjudicate this rank.
            if truth and a is None:
                continue
            self.matrix[rank, a, b] += 1

    def tables(self, sample):
        def div(a, b): return round(a / b, 8) if b else 'NA'
        metrics = []
        confusion = []
        for rank in ('species', 'genus', 'family', 'group'):
            by = defaultdict(lambda: [0, 0, 0])
            micro = [0, 0, 0]
            tested = 0
            for (level, a, b), n in self.matrix.items():
                if level != rank:
                    continue
                tested += n
                confusion.append([sample, rank, a if a is not None else '.', b if b is not None else '.', n])
                if a is not None:
                    by[a][0] += n
                    micro[0] += n
                if b is not None:
                    by[b][1] += n
                    micro[1] += n
                if a is not None and a == b:
                    by[a][2] += n
                    micro[2] += n
            if rank == 'group':
                by['human']
                by['plant']
            for label, (pos, called, tp) in [('ALL', micro), *sorted(by.items(), key=lambda kv: str(kv[0]))]:
                metrics.append([sample, rank, label, self.n, tested, pos, called, tp, called - tp, pos - tp,
                                div(tp, called), div(tp, pos), div(2 * tp, called + pos)])
        return metrics, confusion


def analyze(folder, s, verify, raw_check, truth, single, tax_cache):
    sample = s['sample']
    total = s['pairs']
    bins = {}
    for r in rows(folder / 'taxa.tsv'):
        tid = int(r['taxid'])
        if tid in bins or Path(r['bin']).is_absolute() or '..' in Path(r['bin']).parts:
            raise RuntimeError('Invalid taxon/bin manifest')
        bins[tid] = r
    if sum(int(r['pairs']) for r in bins.values()) != total or len({r['bin'] for r in bins.values()}) != len(bins):
        raise RuntimeError('Taxon bin count/total mismatch: ' + sample)
    for name, size in s['file_sizes'].items():
        if (folder / name).stat().st_size != size:
            raise RuntimeError('Size mismatch: ' + str(folder / name))
    need_assignments = verify or truth
    accuracy = None
    baseline = None
    if truth:
        path = s['taxonomy']
        if path not in tax_cache:
            tax_cache.clear()
            tax_cache[path] = Lineages(Taxonomy(path))
        accuracy = Accuracy(tax_cache[path])
        baseline = Accuracy(tax_cache[path])
    if need_assignments:
        labels = truth_rows(truth, sample, single) if truth else iter(())
        label = next(labels, None)
        seen = Counter()
        reasons = Counter()
        joint = Counter()
        digests = {}
        raw = (fragments(*(x[0] for x in s['request']['raw'])) if raw_check else None)
        whole = new_digest()
        for i, r in enumerate(rows(folder / 'assignments.tsv.gz'), 1):
            if int(r['ordinal']) != i:
                raise RuntimeError('Assignment ordinal mismatch')
            t = int(r['final_taxid'])
            if t not in bins:
                raise RuntimeError('Unknown assignment bin')
            seen[t] += 1
            reasons[r['evidence']] += 1
            joint[t, r['evidence']] += 1
            if raw_check:
                pair = next(raw, None)
                if pair is None:
                    raise RuntimeError('Original FASTQ shorter than assignments')
                if t not in digests:
                    digests[t] = new_digest()
                update_digest(digests[t], *pair)
                update_digest(whole, *pair)
            if label and label[0] == i:
                accuracy.add(label[1], t)
                baseline.add(label[1], int(r['kraken_taxid']))
                label = next(labels, None)
        if seen != {t: int(r['pairs']) for t, r in bins.items()}:
            raise RuntimeError('Assignment totals disagree with bins')
        if reasons != {r['evidence']: int(r['pairs']) for r in rows(folder / 'evidence.tsv')}:
            raise RuntimeError('Assignment evidence totals disagree')
        for r in rows(folder / 'taxon_evidence.tsv'):
            key = int(r['final_taxid']), r['evidence']
            if joint.pop(key, None) != int(r['pairs']):
                raise RuntimeError(f'Assignment taxon/evidence count disagrees: {sample} {key}')
        if joint:
            raise RuntimeError('Missing taxon/evidence rows: ' + sample)
        if label:
            raise ValueError('Truth ordinal outside sample: ' + str(label[0]))
        if truth and not accuracy.n:
            raise ValueError('No matched truth rows for ' + sample)
        if raw_check and (next(raw, None) is not None or whole.hexdigest() != s['input_digest'] or
                          any(digests[t].hexdigest() != r['pairs_sha256'] for t, r in bins.items())):
            raise RuntimeError('Raw reads or assignment-to-bin order disagrees')
    quality = []
    if verify:
        # Phred+33; translate quality bytes to one indicator per base in C.
        q30 = bytes(int(i >= 63) for i in range(256))
        for r in bins.values():
            digest = new_digest()
            n = 0
            totals = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
            base = folder / r['bin']
            for a, b in fragments(base / 'R1.fastq.gz', base / 'R2.fastq.gz'):
                update_digest(digest, a, b)
                for acc, record in zip(totals, (a, b)):
                    seq = record[1].rstrip(b'\r\n').upper()
                    q = record[3].rstrip(b'\r\n')
                    acc[0] += len(seq)
                    acc[1] += seq.count(b'G') + seq.count(b'C')
                    acc[2] += seq.count(b'N')
                    acc[3] += q.translate(q30).count(1)
                    acc[4] += sum(q) - 33 * len(q)
                n += 1
            if n != int(r['pairs']) or digest.hexdigest() != r['pairs_sha256']:
                raise RuntimeError('Physical bin mismatch: ' + str(base))
            for mate, (bases, gc, ns, high, phred) in enumerate(totals, 1):
                quality.append([sample, r['taxid'], r['group'], mate, n, bases,
                               gc / bases, ns / bases, high / bases, phred / bases])
    groups = Counter()
    levels = Counter()
    for r in bins.values():
        n = int(r['pairs'])
        groups[r['group']] += n
        levels[r['rank']] += n
    if sum(int(r['pairs']) for r in rows(folder / 'transitions.tsv')) != total:
        raise RuntimeError('Transition total mismatch')
    return bins, groups, levels, accuracy, baseline, quality


def requested_datasets(params):
    value = params.get('dat_list')
    if value is None:
        return []
    names = [name.strip() for name in value.split(',')]
    for name in names:
        _name(name)
    return list(dict.fromkeys(names))


def discover(cfg, params):
    if params.get('result_dir'):
        p = Path(params['result_dir']).absolute()
        if not (p / 'summary.json').is_file():
            raise ValueError('Result directory needs ranked summary.json')
        return [p]
    root = Path(cfg['output_dir']) / 'ranked'
    names = set(params.get('sample') or [])
    cohorts = requested_datasets(params)
    for name in names:
        parts = name.split('/')
        if len(parts) != 2:
            raise ValueError('Use DATASET/SAMPLE')
        for part in parts:
            _name(part)
        if cohorts and parts[0] not in cohorts:
            raise ValueError('Sample outside selected datasets')
    if names:
        return [(root / name).absolute() for name in sorted(names)]
    datasets = [root / name for name in sorted(cohorts)] if cohorts else sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith('.'))
    result = []
    for ds in datasets:
        if not ds.is_dir():
            if cohorts and (Path(cfg['reads_raw']) / ds.name).is_dir():
                print('No published samples yet: ' + ds.name, flush=True)
                continue
            raise ValueError('Missing ranked cohort: ' + str(ds))
        samples = sorted(p for p in ds.iterdir() if not p.name.startswith('.')
                         and (p / 'summary.json').is_file())
        if not samples:
            print('No published samples yet: ' + ds.name, flush=True)
        result.extend(p.absolute() for p in samples)
    if not result:
        raise ValueError('No completed ranked outputs')
    return result


@contextmanager
def read_sample(path):
    if path.is_symlink():
        raise ValueError('Expected a regular ranked sample directory: '+str(path))
    lock = path.parent / ('.' + path.name + '.runs') / 'publish.lock'
    if lock.is_file():
        from workflow_lock import flock
        with lock.open('rb') as handle:
            flock(handle, fcntl.LOCK_SH)
            yield
    else:
        yield


def index_page(root):
    main = root / 'main'
    stage = main / '.index.tmp'
    stage.write_text(
        '<!doctype html><html lang="zh"><meta charset="utf-8"><title>METAG</title>'
        '<meta http-equiv="refresh" content="0;url=ranked/explore.html">'
        '<a href="ranked/explore.html">打开交互报告</a></html>', encoding='utf-8')
    os.replace(stage, main / 'index.html')


def report(cfg, params):
    root = Path(params['out'])
    main = root / 'main'
    dest = main / 'ranked'
    folders = discover(cfg, params)
    for path in folders:
        if root.resolve().is_relative_to(path.resolve()) or path.resolve().is_relative_to(root.resolve()):
            raise ValueError('Choose a separate paper output root')
    with staging(main, '.building-') as folder:
        write_report(folder, folders, params, params.get('check_bins') or params.get('raw_check'),
                     refine_root=Path(cfg['work_dir'])/'refine')
        script = os.environ.get('METAG_RSCRIPT')
        if script:
            runner = [script]
        elif resolve_mode() == 'hpc':
            runner = ['bash', str(SHARED_RUNTIME / '0os.r.sh')]
        else:
            runner = ['Rscript']
        subprocess.run([*runner, '--vanilla', str(CODE_ROOT / 'f' / 'metag.paper.R'),
                        str(folder), *requested_datasets(params)], check=True)
        seal_audit.summarize_windows(folder)
        result = render(folder)
        explorer.write(folder)
        with locked(root / 'logs/publish.lock'):
            dest.parent.mkdir(parents=True, exist_ok=True)
            recover_directory(dest)
            if dest.is_symlink() or (dest.exists() and not flat_ranked(dest)):
                raise ValueError('Report destination is already in use: ' + str(dest))
            publish_directory(folder, dest)
            index_page(root)
    print(f'{len(result)} sample(s): {dest}/report.html', flush=True)


def write_evidence(path, sample, bins, writers):
    """Check both margins before reporting how each target species bin was assigned."""
    by_reason = Counter()
    by_taxon = Counter()
    for r in rows(path / 'taxon_evidence.tsv'):
        tid, reason, n = int(r['final_taxid']), r['evidence'], int(r['pairs'])
        if n < 0 or tid not in bins:
            raise ValueError(f'Invalid evidence row for {sample}: {tid}, {reason}, {n}')
        by_taxon[tid] += n
        by_reason[reason] += n
        writers['taxon_evidence'].writerow([sample, tid, reason, n])
        taxon = bins[tid]
        if taxon['rank'] == 'species' and taxon['group'] in ('human', 'plant', 'virus'):
            writers['target_species_evidence'].writerow(
                [sample, taxon['group'], tid, taxon['name'], reason, n])
    declared = {}
    for r in rows(path / 'evidence.tsv'):
        reason, n = r['evidence'], int(r['pairs'])
        if reason in declared or n < 0:
            raise ValueError(f'Invalid evidence total for {sample}: {reason}')
        declared[reason] = n
        writers['evidence'].writerow([sample, reason, n])
    if by_taxon != {tid: int(r['pairs']) for tid, r in bins.items()} or by_reason != declared:
        raise ValueError('Taxon evidence disagrees with bin or evidence totals: ' + sample)


def write_report(out, folders, p, verify, refine_root=None):
    from contextlib import ExitStack
    from collections import OrderedDict
    truth = p.get('truth')
    header = 'sample pairs species_pairs higher_rank_pairs unclassified_pairs human_pairs plant_pairs virus_pairs integrity truth_pairs database code rank_min_windows'.split()
    schema = {
        'samples': header, 'taxa': 'sample taxid name group rank pairs bin'.split(),
        'groups': 'sample group pairs total fraction'.split(), 'ranks': 'sample rank pairs total fraction'.split(),
        'evidence': 'sample evidence pairs'.split(), 'transitions': 'sample kraken_group final_group pairs'.split(),
        'taxon_evidence': 'sample final_taxid evidence pairs'.split(),
        'target_species_evidence': 'sample group taxid name evidence pairs'.split(),
        'panels': 'sample taxid name role kraken_pairs near_relative_and_background agreeing_nuclear_windows_capped'.split(),
        'panel_omitted': 'sample taxid kraken_pairs reason'.split()}
    if refine_root is not None:
        schema.update(refine_tables.SCHEMA)
        schema.update(seal_bins=seal_audit.SUMMARY,seal_loci=seal_audit.LOCI,seal_status=["sample","state","panel_id"])
    if truth:
        schema['truth_metrics'] = 'sample rank label labeled evaluable positives called tp fp fn precision recall f1 method'.split()
        schema['truth_confusion'] = 'sample rank truth prediction pairs method'.split()
    if verify:
        schema['fastq_qc'] = 'sample taxid group mate pairs bases gc_fraction n_fraction q30_fraction mean_phred'.split()
    summaries = []
    tax_cache = {}
    snapshots = {}
    for path in folders:
        with read_sample(path):
            s = json.loads((path / 'summary.json').read_text())
        if s.get('format') != FORMAT:
            raise ValueError('Requires ranked output: ' + str(path))
        if not p.get('result_dir') and s['sample'] != path.parent.name + '/' + path.name:
            raise ValueError('Sample receipt disagrees with its directory: ' + str(path))
        if s['sample'] in snapshots:
            raise ValueError('Duplicate sample ' + s['sample'])
        snapshots[s['sample']] = path
    with ExitStack() as stack:
        writers = {}
        for title, fields in schema.items():
            h = stack.enter_context((out / (title + '.tsv')).open('w', encoding='utf-8', newline=''))
            writers[title] = csv.writer(h, delimiter='\t', lineterminator='\n')
            writers[title].writerow(fields)
        truth_paths = {}
        if truth and len(folders) > 1:
            # One scan of a cohort truth table, at most 24 descriptors.
            work = Path(stack.enter_context(tempfile.TemporaryDirectory(dir=out)))
            truth_paths = {name: work / str(i) for i, name in enumerate(snapshots)}
            streams = OrderedDict()
            try:
                for row in rows(truth):
                    name = row.get('sample')
                    if not name:
                        raise ValueError('Multi-sample truth requires a sample column')
                    if name not in truth_paths:
                        continue
                    h = streams.pop(name, None)
                    if h is None:
                        if len(streams) >= 24:
                            streams.popitem(last=False)[1].close()
                        first = not truth_paths[name].exists()
                        h = truth_paths[name].open('a')
                        if first:
                            h.write('ordinal\ttrue_taxid\n')
                    streams[name] = h
                    h.write(f"{int(row['ordinal'])}\t{int(row['true_taxid'])}\n")
            finally:
                for h in streams.values():
                    h.close()
        for name, path in snapshots.items():
            with read_sample(path):
                s = json.loads((path / 'summary.json').read_text())
                if s.get('sample') != name or s.get('format') != FORMAT:
                    raise ValueError('Sample changed during report generation: ' + name)
                if refine_root is not None:
                    split_sha=sha256(path/'summary.json')
                    refine_tables.collect(refine_root,name,split_sha,writers)
                    seal_audit.collect(refine_root,name,split_sha,writers)
                labels = truth_paths.get(name, truth)
                if labels and not Path(labels).exists():
                    raise ValueError('Missing truth for ' + name)
                bins, groups, levels, accuracy, baseline, stats = analyze(
                    path, s, verify, p.get('raw_check'), labels, True, tax_cache)
                total = s['pairs']
                item = [name, total, s['exact_species_pairs'], total - s['exact_species_pairs'] - s['unclassified_pairs'],
                        s['unclassified_pairs'], groups['human'], groups['plant'], groups['virus'],
                        'raw_and_bins' if p.get('raw_check') else 'bins' if verify else 'sizes_and_counts',
                        accuracy.n if accuracy else 'NA', s['request']['database'], s['request']['code'],
                        s['request'].get('settings', {}).get('rank_min_windows', 'NA')]
                writers['samples'].writerow(item)
                summaries.append(item[:10])
                if verify:
                    writers['fastq_qc'].writerows(stats)
                for r in bins.values():
                    writers['taxa'].writerow([name, *[r[k] for k in schema['taxa'][1:]]])
                for dest, counts in (('ranks', levels), ('groups', groups)):
                    writers[dest].writerows([name, k, n, total, n / total if total else 'NA']
                                            for k, n in sorted(counts.items()))
                write_evidence(path, name, bins, writers)
                for title, file in [('transitions', 'transitions'), ('panels', 'panel'), ('panel_omitted', 'panel_omitted')]:
                    writers[title].writerows([name, *[r[k] for k in schema[title][1:]]]
                                             for r in rows(path / (file + '.tsv')))
                if accuracy:
                    for method, scorer in [('full_kraken', baseline), ('final', accuracy)]:
                        metrics, confusion = scorer.tables(name)
                        writers['truth_metrics'].writerows([*r, method] for r in metrics)
                        writers['truth_confusion'].writerows([*r, method] for r in confusion)
    return summaries


def paper_release(cfg):
    return fingerprint([release(dict(cfg))[0], sha256(CODE_ROOT / 'f' / 'metag.paper.R'), refine_tables.code_id(),
                        [(p.name,sha256(p)) for p in sorted((CODE_ROOT/'paper').glob('explorer.*'))]])


def execute_paper(job):
    log = Path(job['log'])
    info = dict(action=job['params']['action'], started=datetime.now(timezone.utc).isoformat())
    with (log/'report.log').open('w', encoding='utf-8') as stream, redirect_stdout(stream), redirect_stderr(stream):
        atomic_json(log/'status.json', dict(info, status='running'))
        try:
            if job.get('code') and paper_release(job['config']) != job['code']:
                raise RuntimeError('Scripts changed after planning. Submit the report again.')
            report(job['config'], job['params'])
        except BaseException:
            traceback.print_exc()
            atomic_json(log/'status.json', dict(info, status='failed', finished=datetime.now(timezone.utc).isoformat()))
            raise
        atomic_json(log/'status.json', dict(info, status='complete', finished=datetime.now(timezone.utc).isoformat()))
    print(Path(job['params']['out'])/'main/index.html', flush=True)


def main():
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        job = json.loads(Path(sys.argv[2]).read_text())
        execute_paper(job)
        return
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('run',))
    parser.add_argument('--sample', action='append')
    for k in ('dat-list', 'result-dir', 'config', 'root', 'truth', 'out', 'queue'):
        parser.add_argument('--' + k)
    parser.add_argument('--check-bins', action='store_true')
    parser.add_argument('--raw-check', action='store_true')
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--submit', choices=('local', 'hpc', 'auto'), default=None)
    for k in ('cores', 'mem-gib', 'hours'):
        parser.add_argument('--' + k, type=positive)
    a = parser.parse_args()
    if a.result_dir and (a.sample or a.dat_list):
        parser.error('--result-dir excludes dataset/sample selectors')
    mode = resolve_mode(a.submit)
    os.environ['submit'] = os.environ['SUBMIT'] = mode
    cfg = load_config(a.config, a.root)
    params = vars(a).copy()
    for k in ('result_dir', 'truth', 'out'):
        if params[k]:
            params[k] = str(Path(params[k]).absolute())
    root = Path(params['out'] or Path(cfg['work_dir'])/'paper').absolute()
    if params['result_dir']:
        inp = Path(params['result_dir']).resolve()
        if root.resolve().is_relative_to(inp) or inp.is_relative_to(root.resolve()):
            parser.error('--out and --result-dir must be separate')
    root = paper_root(root)
    params['out'] = str(root)
    cfg.update(classify_cores=a.cores or 2, classify_memory_gib=a.mem_gib or 32,
               classify_hours=a.hours or 24, classify_queue=a.queue or cfg['classify_queue'])
    # Data discovery still uses the original config. Only the scheduler lives here.
    scheduler_cfg = dict(cfg, work_dir=str(root), code_root=str(CODE_ROOT))
    folder = new_run(root / 'jobs')
    log = root / 'logs' / folder.name
    log.mkdir()
    job = dict(config=cfg, params=params, log=str(log), code=paper_release(cfg))
    spec = folder / 'request.json'
    atomic_json(spec, job)
    cmd = write_cmd(scheduler_cfg, folder/'paper.cmd', 'classify',
                    [sys.executable, '-B', str(CODE_ROOT/'f/rank_paper.py'), '--worker', str(spec)], 'metag.paper')
    print(cmd, flush=True)
    if not a.plan_only:
        if mode == 'local':
            execute_paper(job)
        else:
            submit_cmd(scheduler_cfg, cmd)


if __name__ == '__main__':
    main()
