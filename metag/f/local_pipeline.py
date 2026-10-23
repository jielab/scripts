#!/usr/bin/env python3
"""Resumable species splitting of externally checked FASTQ inputs."""
import csv
import copy
import gzip
import hashlib
import html
import io
from contextlib import nullcontext
from functools import lru_cache
import multiprocessing
import os
import signal
import stat
import time
import json
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from metag import json_write, write_tsv, now, lock, tool_version
from splitter import GzipBins, split_stream, ASSIGNMENT_FIELDS
from taxonomy import Taxonomy, prepare_sqlite, snapshot
from raw_cleanup import cleanup_dataset

VERSION = 2
DATASET_SCHEMA = 'metag-dataset-split-v1'
NAME = re.compile(r'(.+?)(?:_([12])|_R([12])(?:_001)?)?\.(?:fastq|fq)(?:\.gz)?$', re.I)


def identity(paths):
    rows = []
    for p in sorted(paths):
        info = p.stat()
        rows.append(dict(path=str(p.absolute()), target=str(p.resolve()),
                         bytes=info.st_size, mtime_ns=info.st_mtime_ns))
    return rows


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def code_identity():
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).parent.glob('*.py'))}


def split_directory(args):
    return (getattr(args, 'dir_split', None) or args.dir_fastq.resolve().parent / 'split').resolve()


def discover(args):
    raw = args.dir_fastq.resolve()
    if not raw.is_dir():
        raise RuntimeError(f'FASTQ 目录不存在：{raw}')
    finished = set()
    split = split_directory(args)
    if split.is_dir():
        finished = {p.name for p in split.iterdir()
                    if p.is_dir() and not p.is_symlink() and (p / 'split.done').is_file()}
    names = ([x.strip() for x in args.dat_list.split(',')] if args.dat_list is not None
             else sorted(finished | {p.name for p in raw.iterdir()
                                     if p.is_dir() and not p.is_symlink()}))
    if not names or len(names) != len(set(names)):
        raise RuntimeError('没有数据集目录，或 --dat-list 存在重复项')
    for name in names:
        if not name or name in ('.', '..') or '/' in name or '\\' in name:
            raise RuntimeError(f'无效数据集名称：{name!r}')
        if (raw / name).is_symlink() or (not (raw / name).is_dir() and name not in finished):
            raise RuntimeError(f'数据集目录不存在或是目录链接：{raw / name}')
    return raw, names


def entries(root):
    # scandir-backed walk avoids an extra is_dir/stat call for every FASTQ on WSL.
    # Do not traverse directory links, but report them (including broken links).
    found = []
    pending = [root]
    visited = 0
    counts = Counter()
    def progress(final=False):
        print(f'[scan] {"finish" if final else "progress"} {root.name}: '
              f'directories={visited:,}; FASTQ candidates={counts["fastq"]:,}; '
              f'other files={counts["other"]:,}; links={counts["links"]:,} '
              '(file contents are not pre-read)', flush=True)
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as items:
            for item in items:
                if item.is_symlink():
                    counts['links'] += 1
                    found.append(Path(item.path))
                elif item.is_dir(follow_symlinks=False):
                    pending.append(Path(item.path))
                else:
                    counts['fastq' if NAME.fullmatch(item.name) and item.is_file(follow_symlinks=False) else 'other'] += 1
                    found.append(Path(item.path))
        visited += 1
        if visited % 1000 == 0:
            progress()
    progress(final=True)
    return sorted(found)


def inventory(raw, names, manifest=None):
    """Track selected FASTQ inputs, never ancillary files or ignored links.

    During a run the input manifest fixes the inputs. On resume, discover
    current FASTQ candidates as well so newly added runs are included.
    """
    raw = raw.resolve()
    if manifest is None:
        paths = [p for name in names for p in entries(raw / name)
                 if NAME.fullmatch(p.name) and p.is_file()]
    else:
        selected = set(names)
        paths = [Path(row['path']) for row in manifest if row['dataset'] in selected]
    rows = []
    for p in sorted(set(paths)):
        info = p.lstat()
        target, target_info = str(p), info
        if stat.S_ISLNK(info.st_mode):
            target = str(p.resolve())
            try:
                target_info = p.stat()
            except FileNotFoundError:
                target_info = None
        rows.append((str(p.relative_to(raw)), info.st_size, info.st_mtime_ns,
                     target, target_info.st_size if target_info else None,
                     target_info.st_mtime_ns if target_info else None))
    return rows


def report(folder, title, counts):
    """Dependency-free static overview; exact values also published as TSV."""
    folder.mkdir(parents=True, exist_ok=True)
    write_tsv(folder / 'summary.tsv', ['category', 'count'],
              [dict(category=k, count=v) for k, v in counts.items()])
    height = 90 + 36 * len(counts)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="{height}" viewBox="0 0 1000 {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="20" y="30" font-size="20">{html.escape(title)}</text>']
    maximum = max(counts.values(), default=1) or 1
    for i, (key, value) in enumerate(counts.items()):
        y = 65 + i * 36
        svg.extend([f'<text x="20" y="{y+17}">{html.escape(str(key))}</text>',
                    f'<rect x="300" y="{y}" width="{550*value/maximum}" height="24" fill="#367baa"/>',
                    f'<text x="870" y="{y+17}">{value:,}</text>'])
    (folder / 'summary.svg').write_text('\n'.join(svg + ['</svg>']))


def outputs(root):
    return identity(p for p in root.rglob('*') if p.is_file() and p.name != 'split.done')


def worker_start():
    # A process group per worker also contains its Kraken2/decompression children.
    os.setsid()
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def parallel_results(function, tasks, jobs):
    """Spawn avoids inheriting workflow locks and SQLite connections."""
    if jobs == 1:
        for task in tasks:
            yield function(task)
        return
    pool = multiprocessing.get_context('spawn').Pool(jobs, initializer=worker_start)
    try:
        yield from pool.imap_unordered(function, tasks, chunksize=1)
        pool.close()
        pool.join()
    except BaseException:
        # Terminate process groups, including subprocesses, before joining the pool.
        groups = [worker.pid for worker in pool._pool]
        # Let Pool stop its queue-management threads before killing workers;
        # external kills first can leave its queue lock held and deadlock terminate.
        pool.terminate()
        for pgid in groups:
            try:
                os.killpg(pgid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        pool.join()
        raise


# 🚩 Input discovery without reading FASTQ contents
def input_files(raw, names, previous=()):
    """Trust external MD5 checks; only inspect paths, sizes and filename roles."""
    manifest = []
    for dataset in names:
        print(f'[metag] START inputs dataset={dataset}; external MD5 assumed', flush=True)
        groups = defaultdict(list)
        ignored = 0
        for p in entries(raw / dataset):
            if (p.name.endswith(('.aria2', '.tmp'))
                    and NAME.fullmatch(p.name.rsplit('.', 1)[0])):
                raise RuntimeError(f'FASTQ 下载尚未完成，保留原文件：{p}')
            match = NAME.fullmatch(p.name)
            if not match:
                ignored += 1
                continue
            try:
                info = p.stat()
            except FileNotFoundError as error:
                raise RuntimeError(f'FASTQ 缺失或链接失效：{p}') from error
            if not stat.S_ISREG(info.st_mode):
                raise RuntimeError(f'FASTQ 缺失或不是普通文件（含失效链接）：{p}')
            if info.st_size == 0:
                raise RuntimeError(f'空 FASTQ 文件：{p}')
            run = str(p.relative_to(raw / dataset).parent / match[1])
            role = 'r' + (match[2] or match[3]) if match[2] or match[3] else 'single'
            groups[run].append(dict(dataset=dataset, run=run, role=role, path=str(p)))
        if not groups:
            raise RuntimeError(f'数据集没有可用 FASTQ：{raw / dataset}')
        for run, rows in sorted(groups.items()):
            roles = [row['role'] for row in rows]
            paths = ', '.join(row['path'] for row in rows)
            if len(roles) != len(set(roles)):
                raise RuntimeError(f'FASTQ 角色重复、无法确定配对：{dataset}/{run}；{paths}')
            if ('r1' in roles) != ('r2' in roles):
                raise RuntimeError(f'FASTQ 缺失 mate：{dataset}/{run}；现有文件：{paths}')
            manifest.extend(rows)
        print(f'[metag] DONE inputs dataset={dataset}; files={sum(map(len, groups.values()))}; '
              f'runs={len(groups)}; ignored={ignored}; no content pre-scan', flush=True)
    present = {row['path'] for row in manifest}
    # Retain knowledge of whole runs that disappear between submissions. Finished
    # datasets are excluded by the caller and may have had their raw removed.
    for row in previous:
        dataset = row['dataset']
        p = Path(row['path'])
        if (dataset in names and p.is_relative_to(raw / dataset)
                and NAME.fullmatch(p.name) and str(p) not in present):
            raise RuntimeError(f'上次输入清单中的 FASTQ 缺失：{p}')
    return sorted(manifest, key=lambda row: (row['dataset'], row['run'], row['role']))


@lru_cache(maxsize=50000)
def bin_directory(group, bucket):
    prefix = 'species' if re.match(r'^[1-9][0-9]*_', bucket) else 'split.bad'
    leaf = bucket.rsplit('/', 1)[-1]
    taxid = '0' if leaf == 'taxid_0' else leaf.split('_', 1)[0]
    return f'{prefix}/{group}/{taxid}'


class RoutedBins:
    def __init__(self, root, tax):
        self.writer = GzipBins(root)
        self.tax = tax

    def write(self, key, record):
        # Taxonomy's species bucket starts with a numeric taxid; unresolved buckets do not.
        group, rest = key.split('/', 1)
        bucket, filename = rest.rsplit('/', 1)
        self.writer.write(f'{bin_directory(group, bucket)}/{filename}', record)

    def close(self):
        self.writer.close()


def classify_run(args, rows, root, tax, db, signature):
    done = root / 'split.done'
    if done.is_file():
        state = json.loads(done.read_text())
        if state['signature'] != signature:
            raise RuntimeError(f'输入/参数/参考发生变化；请使用新的 --dir-split，保留旧结果：{root}')
        if state['outputs'] != outputs(root):
            raise RuntimeError(f'已完成输出缺失或发生变化：{root}')
        print(f'[split] skip {root.name}', flush=True)
        return state['counts']
    if root.exists():
        raise RuntimeError(f'发现没有完成标记的正式输出：{root}')
    root.parent.mkdir(parents=True, exist_ok=True)
    stage = root.with_name(root.name + '.partial')
    owner = dict(kind='metag-split-stage-v2', final_root=str(root))
    if stage.exists() or stage.is_symlink():
        receipt = stage / '.metag-stage.json'
        if stage.is_symlink() or not receipt.is_file() or json.loads(receipt.read_text()) != owner:
            raise RuntimeError(f'拒绝清理无法确认来源的暂存目录：{stage}')
        # Only our own unfinished output is disposable; never remove final results or inputs.
        shutil.rmtree(stage)
    stage.mkdir()
    json_write(stage / '.metag-stage.json', owner)
    bins = RoutedBins(stage, tax)
    counter = Counter()
    paths = {r['role']: Path(r['path']) for r in rows}
    input_state = identity(paths.values())
    jobs = []
    if 'r1' in paths:
        jobs.append(('PE', paths['r1'], paths['r2']))
    if 'single' in paths:
        jobs.append(('SINGLE', paths['single'], None))
    try:
        for stream, r1, r2 in jobs:
            cmd = ['kraken2', '--db', str(db), '--threads', str(args.threads),
                   '--confidence', str(args.confidence), '--minimum-hit-groups', str(args.hit_groups),
                   '--report', str(stage / f'{stream}.kreport'), '--report-minimizer-data']
            if args.classification_mode == 'quick':
                cmd.append('--quick')
            if args.memory_mapping:
                cmd.append('--memory-mapping')
            if r2:
                cmd += ['--paired', str(r1), str(r2)]
            else:
                cmd.append(str(r1))
            json_write(stage / f'{stream}.command.json', cmd)
            assignment_file = (io.BufferedWriter(gzip.open(stage / f'{stream}.assignments.tsv.gz',
                               'wb', compresslevel=1), buffer_size=1 << 20)
                               if args.keep_assignments else nullcontext(None))
            with (stage / f'{stream}.log').open('w') as err, assignment_file as assignments:
                if assignments is not None:
                    assignments.write(('\t'.join(ASSIGNMENT_FIELDS) + '\n').encode())
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err)
                try:
                    c, totals = split_stream(proc.stdout, r1, r2, stream, tax, bins, assignments)
                    if proc.wait():
                        raise RuntimeError(f'Kraken2 失败：{stage / (stream + ".log")}')
                    counter.update(c)
                except BaseException:
                    if proc.poll() is None:
                        proc.terminate()
                        try:
                            proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            proc.kill()
                            proc.wait()
                    raise
                finally:
                    proc.stdout.close()
        bins.close()
        if input_state != identity(paths.values()):
            raise RuntimeError('拆分期间输入文件发生变化')
        tax.virus.assert_unchanged()
        counts = [dict(group=g, bucket=b, species_taxid=s, species_name=n, stream=t,
                       fragments=count, reads=count * mates,
                       destination='species' if s else 'split.bad', relative_dir=bin_directory(g, b),
                       reason='' if s else ('kraken_unclassified' if b == 'taxid_0' else 'no_species_in_taxonomic_lineage'))
                  for (g, b, s, n, t, mates), count in sorted(counter.items())]
        write_tsv(stage / 'counts.tsv', ['group', 'bucket', 'species_taxid', 'species_name', 'stream', 'fragments', 'reads', 'destination', 'relative_dir', 'reason'], counts)
        # Always log unresolved output, independently of --keep-assignments.
        # Each bin has one routing reason, applying to all records in its files.
        with (stage / 'split.bad.jsonl').open('w') as log:
            for row in counts:
                if row['destination'] != 'split.bad':
                    continue
                stream = row['stream']
                filenames = ([f'{stream}.R1.fastq.gz', f'{stream}.R2.fastq.gz']
                             if stream == 'PE' else [f'{stream}.fastq.gz'])
                record = dict(row, time=now(),
                              input_files=[str(paths[r]) for r in (('r1', 'r2') if stream == 'PE' else ('single',))],
                              output_files=[str(root / row['relative_dir'] / name) for name in filenames])
                log.write(json.dumps(record, ensure_ascii=False) + '\n')
                print(f'[split.bad] {root.name}/{row["relative_dir"]}/{stream}: '
                      f'{row["reads"]} reads; reason={row["reason"]}', flush=True)
        json_write(stage / 'input.json', dict(files=rows, signature=signature))
        # Record final paths while still staged; publish marker and files atomically together.
        final_outputs = outputs(stage)
        for item in final_outputs:
            for key in ('path', 'target'):
                item[key] = str(root / Path(item[key]).relative_to(stage))
        json_write(stage / 'split.done', dict(signature=signature, completed=now(), counts=counts,
                                            outputs=final_outputs))
        stage.rename(root)
        return counts
    except BaseException as error:
        bins.close()
        failure = f'拆分失败（输入：{", ".join(map(str, paths.values()))}）：{error}'
        if stage.exists():
            json_write(stage / 'FAILED.json', dict(error=failure, time=now()))
        if isinstance(error, (OSError, RuntimeError, ValueError)):
            raise RuntimeError(failure) from error
        raise


def run_stem(runid):
    # A flat accession is already unique within a dataset; nested/general names
    # receive a short stable suffix without copying their directory ancestry.
    if re.fullmatch(r'[SED]RR[0-9]+', runid):
        return runid
    return re.sub(r'[^A-Za-z0-9_.-]', '_', Path(runid).name)[:40] + '-' + digest(runid)[:12]


def classify_task(task):
    args, rows, root, taxonomy_path, mapping, db, signature, dataset, runid = task
    print(f'[split] start {dataset}/{runid}', flush=True)
    tax = Taxonomy(taxonomy_path, mapping)
    try:
        return dataset, runid, classify_run(args, rows, root, tax, db, signature)
    finally:
        tax.close()


def write_split_summary(folder, rows):
    rows.sort(key=lambda r: (r['dataset'], r['run'], r['group'], r['bucket'], r['stream']))
    write_tsv(folder / 'species_by_run.tsv',
              ['dataset', 'run', 'group', 'bucket', 'species_taxid', 'species_name',
               'stream', 'fragments', 'reads', 'destination', 'relative_dir', 'reason'], rows)
    totals = Counter()
    for row in rows:
        totals[row['destination'] + '/' + row['group']] += row['reads']
    report(folder, 'Kraken2 read assignments (not abundance)', totals)
    write_tsv(folder / 'virus_by_run.tsv',
              ['dataset', 'run', 'group', 'species_taxid', 'species_name', 'reads',
               'destination', 'relative_dir'],
              [r for r in rows if r['group'].startswith('virus_')])


def dataset_receipt(root, raw, dataset):
    marker = root / 'split.done'
    if not marker.is_file():
        return None
    state = json.loads(marker.read_text())
    if (state.get('schema') != DATASET_SCHEMA or state.get('dataset') != dataset
            or state.get('raw_directory') != str(raw / dataset)
            or not state.get('runs') or not state.get('artifacts')
            or not state.get('inventory') or not isinstance(state.get('context'), dict)):
        raise RuntimeError(f'数据集完成标记不匹配，保留现场：{marker}')
    return state


def restore_dataset(root, raw, dataset, state, context):
    """Validate saved products; a completed dataset no longer needs its raw."""
    if state['context'] != context:
        raise RuntimeError(f'已完成数据集的代码/参数/参考发生变化，请使用新的 --dir-split：{root}')
    expected_paths = [root / item['stem'] / 'split.done' for item in state['runs']]
    expected_paths += sorted((root / 'summary').glob('*'))
    if (any(item['stem'] != run_stem(item['run']) for item in state['runs'])
            or len({item['run'] for item in state['runs']}) != len(state['runs'])
            or state['artifacts'] != identity(expected_paths)):
        raise RuntimeError(f'已完成数据集的完成标记/汇总发生变化：{root}')
    summary = []
    for item in state['runs']:
        run_root = root / item['stem']
        saved = json.loads((run_root / 'split.done').read_text())
        if saved['signature'] != item['signature'] or saved['outputs'] != outputs(run_root):
            raise RuntimeError(f'已完成输出缺失或发生变化：{run_root}')
        summary.extend(dict(dataset=dataset, run=item['run'], **row) for row in saved['counts'])
    # Missing originals are allowed; changed/new FASTQ inputs are not.
    source = raw / dataset
    if source.exists():
        if source.is_symlink() or not source.is_dir():
            raise RuntimeError(f'已完成数据集的原始目录类型发生变化：{source}')
        previous = {row[0]: digest(row) for row in state['inventory']}
        if any(previous.get(row[0]) != digest(row) for row in inventory(raw, [dataset])):
            raise RuntimeError(f'已完成数据集包含新增或修改的 FASTQ，请使用新的 --dir-split：{source}')
    return summary


def announce_dataset(dataset, raw, root, runs, rows, reused=False):
    print(f'[metag] {"SKIP" if reused else "DONE"} split dataset={dataset}; '
          f'runs={runs}; reads={sum(int(r["reads"]) for r in rows):,}; '
          f'output={root}; 已完成校验，原始 FASTQ 已可释放：{raw / dataset}', flush=True)


def check_split_space(split, tasks):
    # A lower-bound planning check, not a promise about gzip compression ratios.
    # Finished runs need no additional FASTQ space when resuming a dataset.
    required = sum(Path(row['path']).stat().st_size
                   for task in tasks if not (task[2] / 'split.done').is_file()
                   for row in task[1])
    free = shutil.disk_usage(split).free
    if required > free:
        raise RuntimeError(f'split 可用空间不足：{split} 剩余 {free / 2**30:.1f} GiB，'
                           f'下一数据集未完成输入为 {required / 2**30:.1f} GiB（输出可能更大）。'
                           '请先移走已显示 DONE 的数据集 raw，再重跑原命令；已完成结果会复用。')


def run(args, prepare_reference):
    raw, names = discover(args)
    out = args.out.resolve()
    split = split_directory(args)
    delete_raw = getattr(args, 'delete_raw_after_split', False)
    if not isinstance(delete_raw, bool):
        raise ValueError('delete_raw_after_split 必须为布尔值 TRUE/FALSE')
    for target in (out, split):
        if target == raw or target.is_relative_to(raw) or raw.is_relative_to(target):
            raise RuntimeError('输出目录必须与 raw 目录分离')
    print(f'数据集：{", ".join(names)}\nFASTQ：{raw}\n输出：{out}\n拆分输出：{split}', flush=True)
    print(f'[metag] delete_raw_after_split={str(delete_raw).upper()}; '
          '按数据集完成标记和输出校验决定是否清理原始 FASTQ', flush=True)
    if args.dry_run:
        print('[inputs] 依赖外部 MD5；仅发现文件并核对路径/配对，不预读 FASTQ 内容')
        prepare_reference()
        print(f'[split] workers={min(args.split_jobs, args.threads)}; total Kraken2 threads={args.threads}; mode={args.classification_mode}; assignments={args.keep_assignments}')
        print(f'[split] 按数据集顺序 -> {split}/<dataset>/<run>/species 或 split.bad；'
              f'逐数据集 split.done 校验后{"自动删除已处理 raw FASTQ" if delete_raw else "保留 raw FASTQ"}；已完成结果校验后复用。')
        return
    out.mkdir(parents=True, exist_ok=True)
    split.mkdir(parents=True, exist_ok=True)
    # Output locks serialize splitting, completion validation and raw cleanup.
    with (lock(out / '.metag.workflow.lock'), lock(split / '.metag.split.lock')):
        saved_datasets = {name: state for name in names
                          if (state := dataset_receipt(split / name, raw, name)) is not None}
        active_names = [name for name in names if name not in saved_datasets]
        inputs = out / 'inputs'
        inputs.mkdir(exist_ok=True)
        marker = inputs / 'ready.done'
        manifest_path = inputs / 'files.tsv'
        # The old QC table is only a list of expected paths, never evidence that
        # contents were checked in this invocation. Leave historical QC intact.
        previous_path = manifest_path if manifest_path.is_file() else out / 'qc/files.tsv'
        previous = []
        if previous_path.is_file():
            with previous_path.open() as f:
                previous = list(csv.DictReader(f, delimiter='\t'))
        marker.unlink(missing_ok=True)
        manifest = input_files(raw, active_names, previous)
        write_tsv(manifest_path, ['dataset', 'run', 'role', 'path'], manifest)
        current = inventory(raw, active_names, manifest)
        # Reference preparation/audit is idempotent and independent of raw FASTQ.
        prepare_reference()
        db = args.dir_annot.resolve() / '0.metadata/kraken2_broad_v1'
        mapping = args.dir_annot.resolve() / '0.metadata/virus_map_v1/virus_genome_types.sqlite'
        reference = snapshot(db)
        reference['mapping'] = identity([mapping])
        if inventory(raw, active_names, manifest) != current:
            raise RuntimeError('参考准备期间参与分析的 FASTQ 发生变化，请核对后重跑原命令')
        json_write(marker, dict(version=VERSION, content_validation='external_md5_assumed',
                                raw=str(raw), datasets=active_names, inventory=current,
                                manifest_sha=digest(manifest_path.read_text()), reference=reference, completed=now()))
        args = copy.copy(args)
        args.split_jobs = min(args.split_jobs, args.threads)
        args.threads = max(1, args.threads // args.split_jobs)
        args.memory_mapping = args.memory_mapping or args.split_jobs > 1
        version = tool_version(['kraken2', '--version'])
        if args.threads > 1 and not re.search(r'(?<![0-9])2\.17\.1(?![0-9.])', version):
            raise RuntimeError('多线程拆分要求 Kraken2 2.17.1；其他版本请用 --threads 1')
        taxonomy_path = prepare_sqlite(db)
        grouped = defaultdict(list)
        for row in manifest:
            grouped[(row['dataset'], row['run'])].append(row)
        summary = []
        (split / 'split.done').unlink(missing_ok=True)
        code = code_identity()
        context = dict(code=code, reference=reference, version=version, confidence=args.confidence,
                       hit_groups=args.hit_groups, keep_assignments=args.keep_assignments,
                       classification_mode=args.classification_mode)
        print(f'[split] {len(grouped)} runs; workers={args.split_jobs}; '
              f'Kraken2 threads/worker={args.threads}; mmap={args.memory_mapping}; mode={args.classification_mode}', flush=True)
        run_count = 0
        for dataset in names:
            root = split / dataset
            if dataset in saved_datasets:
                saved = saved_datasets[dataset]
                rows = restore_dataset(root, raw, dataset, saved, context)
                summary.extend(rows)
                run_count += len(saved['runs'])
                if delete_raw:
                    cleanup_dataset(root, raw, dataset, context)
                announce_dataset(dataset, raw, root, len(saved['runs']), rows, reused=True)
                continue
            expected = [row for row in current if Path(row[0]).parts[0] == dataset]
            if inventory(raw, [dataset], manifest) != expected:
                raise RuntimeError(f'待拆分数据集参与分析的 FASTQ 发生变化：{dataset}')
            tasks = []
            runs = []
            for (name, runid), rows in sorted(grouped.items()):
                if name != dataset:
                    continue
                stem = run_stem(runid)
                sig = digest(dict(context, files=identity(Path(r['path']) for r in rows)))
                tasks.append((args, rows, root / stem, taxonomy_path, mapping, db, sig, dataset, runid))
                runs.append(dict(run=runid, stem=stem, signature=sig))
            check_split_space(split, tasks)
            print(f'[metag] START split dataset={dataset}; runs={len(tasks)}; output={root}', flush=True)
            rows = []
            results = parallel_results(classify_task, tasks, args.split_jobs)
            try:
                for completed, (_, runid, counts) in enumerate(results, 1):
                    rows.extend(dict(dataset=dataset, run=runid, **c) for c in counts)
                    print(f'[split] dataset={dataset} completed {completed}/{len(tasks)}: {runid}', flush=True)
            finally:
                results.close()
            if inventory(raw, [dataset], manifest) != expected:
                raise RuntimeError(f'拆分期间参与分析的 FASTQ 发生变化：{dataset}')
            if (snapshot(db) != {k: v for k, v in reference.items() if k != 'mapping'}
                    or identity([mapping]) != reference['mapping']):
                raise RuntimeError('拆分期间参考库发生变化')
            # Verify all published runs before announcing that originals may go.
            for item in runs:
                run_root = root / item['stem']
                saved = json.loads((run_root / 'split.done').read_text())
                if saved['signature'] != item['signature'] or saved['outputs'] != outputs(run_root):
                    raise RuntimeError(f'已完成输出缺失或发生变化：{run_root}')
            write_split_summary(root / 'summary', rows)
            artifacts = identity([root / item['stem'] / 'split.done' for item in runs]
                                 + list((root / 'summary').iterdir()))
            json_write(root / 'split.done', dict(
                schema=DATASET_SCHEMA, dataset=dataset, raw_directory=str(raw / dataset),
                context=context, inventory=expected, runs=runs, artifacts=artifacts, completed=now()))
            summary.extend(rows)
            run_count += len(runs)
            # Cleanup reads only the completed dataset inventory; later datasets
            # and the final summaries no longer depend on these raw FASTQ files.
            if delete_raw:
                cleanup_dataset(root, raw, dataset, context)
            announce_dataset(dataset, raw, root, len(runs), rows)
        if snapshot(db) != {k: v for k, v in reference.items() if k != 'mapping'}:
            raise RuntimeError('拆分期间参考库发生变化')
        write_split_summary(out / 'summary', summary)
        json_write(split / 'split.done', dict(completed=now(), datasets=names, runs=run_count, reference=reference))


