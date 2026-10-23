#!/usr/bin/env python3
"""Explicit dataset-list operations; each dataset keeps its own metadata and results."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import uuid

from metag import FASTQ_NAME, parse_ena_report, read_tsv, write_tsv, manifest

HERE = Path(__file__).resolve().parent
DATA_ROOT = Path('/mnt/f/metag')
REFERENCE = Path('/mnt/f/annot/kraken2_ictv')


def registry():
    with (HERE / 'datasets.example.tsv').open(encoding='utf-8', newline='') as stream:
        return {row['dataset']: row for row in csv.DictReader(stream, delimiter='\t')}


def selected_datasets(value, known):
    names = [name.strip() for name in value.split(',')]
    if not all(names) or len(names) != len(set(names)):
        raise ValueError('数据集列表不能包含空名称或重复名称。')
    unknown = [name for name in names if name not in known]
    if unknown:
        raise ValueError(f'未知数据集：{", ".join(unknown)}。运行 ./metag.sh datasets 查看目录。')
    return names


def local_files(root):
    # Keep directory-entry type information: per-file stat calls are costly on WSL drives.
    with os.scandir(root) as entries:
        paths = [Path(entry.path) for entry in entries
                 if entry.is_file(follow_symlinks=False) or entry.is_symlink()]
    if (root / 'raw').is_dir():
        pending = [root / 'raw']
        while pending:
            with os.scandir(pending.pop()) as entries:
                for entry in entries:
                    # os.walk checks symlink targets to identify directories, which
                    # triggers thousands of slow stats on Windows-mounted drives.
                    if entry.is_dir(follow_symlinks=False):
                        pending.append(Path(entry.path))
                    else:
                        paths.append(Path(entry.path))
    fastqs = [p for p in paths if p.name.endswith(('.fastq.gz', '.fq.gz'))]
    markers = [p for p in paths if p.name.endswith(('.fastq.gz.aria2', '.fq.gz.aria2'))]
    other = [p for p in fastqs if not FASTQ_NAME.fullmatch(p.name)]
    return fastqs, markers, other


def show_datasets(args, known):
    names = selected_datasets(args.datasets, known) if args.datasets else list(known)
    print('数据集\tFASTQ路径数\t下载标记\t非ENA文件名\t已有清单\t目录')
    for name in names:
        root = args.data_root / name
        if not root.is_dir():
            if args.datasets:
                print(f'{name}\t-\t-\t-\t-\t{root}（不存在）')
            continue
        files, markers, other = local_files(root)
        exists = 'yes' if (root / 'metadata/samples.tsv').is_file() else 'no'
        print(f'{name}\t{len(files)}\t{len(markers)}\t{len(other)}\t{exists}\t{root}')
    print('仅列出登记的数据集；0metadata、url 不作为数据集；参考资料独立保存在 /mnt/f/annot/kraken2_ictv。')
    print('文件计数含 raw/ 链接路径，不代表独立样本数或 MD5 校验通过。')


def command(argv, execute):
    print(shlex.join(str(arg) for arg in argv), flush=True)
    if execute:
        subprocess.run([str(arg) for arg in argv], check=True)


def pipeline(*argv):
    return [sys.executable, HERE / 'metag.py', *argv]


def check_downloads(root):
    files, markers, other = local_files(root)
    if markers:
        examples = ', '.join(p.name for p in markers[:3])
        raise RuntimeError(f'仍有 {len(markers)} 个 .aria2 标记：{examples}；先确认下载完成，不自动删除标记。')
    if other:
        examples = ', '.join(p.name for p in other[:3])
        raise RuntimeError(f'{len(other)} 个 FASTQ 使用非 ENA 文件名（如 {examples}）；'
                           '需先核对 submitted 文件与 run、处理/合并及配对关系，不能按文件名猜测。')
    if not files:
        raise RuntimeError(f'{root} 和 raw/ 下没有 FASTQ 文件。')


def sync_manifest(args, name, root, record):
    """Reinventory in staging; retain annotations and curated sample identities."""
    meta = root / 'metadata'
    if not args.execute:
        print('计划：重新核对本地文件、加入新增 run，备份并保留已有人工字段和 include。')
        return
    meta.mkdir(exist_ok=True)
    cached = [meta / f'{project}.ena.tsv' for project in record['projects'].split(',')]
    use_cache = all(path.is_file() for path in cached)
    if use_cache:
        try:
            for path in cached:
                parse_ena_report(path.read_text(encoding='utf-8-sig'), str(path))
        except (RuntimeError, UnicodeError):
            use_cache = False
    with tempfile.TemporaryDirectory(prefix='.inventory-', dir=root) as folder:
        stage = Path(folder)
        try:
            manifest(SimpleNamespace(root=root, metadata_dir=stage, refresh=False,
                                     dataset=name, ena_tsv=cached if use_cache else None,
                                     project=record['projects'].split(',')))
        except RuntimeError as error:
            audit = stage / 'inventory.audit.json'
            if audit.exists():
                diagnostic = meta / 'inventory.sync-failed.json'
                shutil.copy2(audit, diagnostic)
                raise RuntimeError(f'{error} 同步未发布，检查记录：{diagnostic}') from error
            raise
        fresh = read_tsv(stage / 'samples.tsv')
        old = read_tsv(meta / 'samples.tsv') if (meta / 'samples.tsv').exists() else []
        previous = {row['run_accession']: row for row in old}
        if len(previous) != len(old):
            raise RuntimeError('旧清单存在重复 run，不能自动同步。')
        missing = set(previous) - {row['run_accession'] for row in fresh}
        if missing:
            raise RuntimeError(f'旧清单中的 run 已不在本地，保留原清单，请检查：{sorted(missing)[:5]}')
        merged = []
        for row in fresh:
            prior = previous.get(row['run_accession'])
            if prior:
                # Layout changes may invalidate curated assay/QC interpretation.
                if prior['observed_layout'] != row['observed_layout']:
                    raise RuntimeError(f"{row['run_accession']} 文件布局变化，请审核后显式 prepare --refresh。")
                inventory = {key: row[key] for key in ('valid', 'observed_layout', 'r1', 'r2', 'single')}
                row.update(prior)
                row.update(inventory)
                for role in ('r1', 'r2', 'single'):
                    if row[role]:
                        row[role] = str(root / 'raw' / row['sample_key'] / row['run_accession'] / Path(row[role]).name)
            merged.append(row)
        by_run = {row['run_accession']: row for row in merged}
        files = read_tsv(stage / 'files.tsv')
        for row in files:
            row['raw_path'] = by_run[row['run_accession']][row['role']]
            if '..' in Path(row['raw_path']).parts or not Path(row['raw_path']).parent.resolve().is_relative_to((root / 'raw').resolve()):
                raise RuntimeError(f"Unsafe raw destination: {row['raw_path']}")
        fields = list(dict.fromkeys(key for row in merged for key in row))
        write_tsv(stage / 'samples.tsv', fields, merged)
        write_tsv(stage / 'files.tsv', list(files[0]), files)
        changed = any(not (meta / path.name).exists() or (meta / path.name).read_bytes() != path.read_bytes()
                      for path in (stage / 'samples.tsv', stage / 'files.tsv'))
        if changed and old:
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '.' + uuid.uuid4().hex[:8]
            backup = root / f'metadata.before-sync.{stamp}'
            shutil.copytree(meta, backup)
            print(f'清单备份：{backup}', flush=True)
        for path in stage.iterdir():
            if changed or path.name not in ('samples.tsv', 'files.tsv'):
                os.replace(path, meta / path.name)
        print(f'同步完成：{len(merged)} 个 run；新增 {len(merged) - len(old)} 个；保留人工注释和 include。')


def prepare(args, name, root, record):
    meta = root / 'metadata'
    if args.execute:
        check_downloads(root)
    if getattr(args, 'sync', False):
        sync_manifest(args, name, root, record)
    elif (meta / 'samples.tsv').exists() and not args.refresh:
        print('复用已有 samples.tsv；需要更新下载状态时使用 prepare --refresh（先备份再重建）。')
    else:
        if args.refresh and meta.exists():
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '.' + uuid.uuid4().hex[:8]
            backup = root / f'metadata.before-refresh.{stamp}'
            if args.execute:
                shutil.copytree(meta, backup)
            print(f'{"备份" if args.execute else "计划备份"}：{meta} -> {backup}')
        argv = pipeline('manifest', '--root', root, '--dataset', name)
        projects = record['projects'].split(',')
        cached = [meta / f'{project}.ena.tsv' for project in projects]
        use_cache = all(path.is_file() for path in cached)
        if use_cache:
            try:
                for path in cached:
                    parse_ena_report(path.read_text(encoding='utf-8-sig'), str(path))
            except (RuntimeError, UnicodeError):
                use_cache = False
        if use_cache:
            for path in cached:
                argv.extend(['--ena-tsv', path])
        else:
            for project in projects:
                argv.extend(['--project', project])
        if args.refresh:
            argv.append('--refresh')
        command(argv, args.execute)
    command(pipeline('organize', '--root', root, '--mode', 'symlink'), args.execute)


def pilot(args, root):
    if not args.execute and not (root / 'metadata/samples.tsv').exists():
        verification = '重新执行 FASTQ MD5 校验' if args.check_md5 else '跳过已外部完成的 FASTQ MD5 校验'
        print(f'计划：先 prepare；再选择最多 {args.limit} 个有效 run，{verification}，'
              f'每输入流最多 {args.max_records} records 试跑，结果写入 {root / "pilot_virus_v2"}。')
        return
    rows = read_tsv(root / 'metadata/samples.tsv')
    eligible = [row for row in rows if row['valid'] == 'yes'
                and row['include'].lower() in ('yes', 'true', '1')
                and row['assay_class'] in ('shotgun_DNA', 'RNA_seq')]
    runs = ','.join(row['run_accession'] for row in eligible[:args.limit])
    if not runs:
        raise RuntimeError('没有可试跑的有效 run；先 prepare 并检查 samples.tsv。')
    out = root / 'pilot_virus_v2'
    command(pipeline('db-audit', '--db', args.db, '--virus-map', args.virus_map,
                     '--threads', args.threads), args.execute)
    if args.check_md5:
        command(pipeline('verify', '--root', root, '--runs', runs, '--md5'), args.execute)
    else:
        print('跳过原始 FASTQ 的重复 MD5 检查（由用户在流程外完成）；需要重查时加 --check-md5。', flush=True)
    command(pipeline('run', '--root', root, '--db', args.db, '--virus-map', args.virus_map,
                     '--runs', runs, '--max-records', args.max_records, '--out', out,
                     '--threads', args.threads), args.execute)
    command(pipeline('summary', '--out', out), args.execute)


def full(args, root):
    argv = ['bash', HERE / 'run_full.sh', '--root', root, '--db', args.db,
            '--virus-map', args.virus_map, '--threads', args.threads]
    if args.execute:
        argv.append('--execute')
    if args.check_md5:
        argv.append('--check-md5')
    command(argv, args.execute)


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('必须为正整数')
    return number


def parser():
    p = argparse.ArgumentParser(prog='./metag.sh', description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    for name, desc in [('datasets', '只读查看目录与下载标记'),
                       ('prepare', '建立/复用清单并整理 raw 软链接'),
                       ('pilot', '选择少量 run 试跑；默认不重复检查 FASTQ MD5'),
                       ('full', '依次全量处理显式选择的数据集')]:
        a = sub.add_parser(name, help=desc, description=desc)
        a.add_argument('datasets', **({'nargs': '?'} if name == 'datasets' else {}),
                       help='英文逗号分隔目录名，例如 emp,Weurope；也支持中文目录名')
        a.add_argument('--data-root', type=Path, default=DATA_ROOT, help='默认 /mnt/f/metag')
        if name == 'datasets':
            continue
        a.add_argument('--execute', action='store_true', help='真正执行；默认只打印计划')
        if name == 'prepare':
            mode = a.add_mutually_exclusive_group()
            mode.add_argument('--sync', action='store_true', help='同步本地文件与新增 run，备份并保留人工注释')
            mode.add_argument('--refresh', action='store_true',
                           help='先备份再重建清单；会重置人工注释，不用于已人工整理的清单')
        else:
            a.add_argument('--db', type=Path, default=REFERENCE / '0.metadata/kraken2_broad_v1',
                           help='默认 /mnt/f/annot/kraken2_ictv/0.metadata/kraken2_broad_v1')
            a.add_argument('--virus-map', type=Path, default=REFERENCE / '0.metadata/virus_map_v1/virus_genome_types.sqlite',
                           help='默认 /mnt/f/annot/kraken2_ictv/0.metadata/virus_map_v1/virus_genome_types.sqlite')
            a.add_argument('--threads', type=positive, default=16)
            a.add_argument('--check-md5', action='store_true',
                           help='重新完整检查原始 FASTQ 的 ENA MD5；默认依赖用户已有的外部校验，不重复读取')
        if name == 'pilot':
            a.add_argument('--limit', type=positive, default=2, help='每数据集最多选择的 run 数，默认 2')
            a.add_argument('--max-records', type=positive, default=100000, help='每输入流上限，默认 100000')
    return p


def main():
    args = parser().parse_args()
    try:
        known = registry()
        args.data_root = args.data_root.resolve()
        if args.command == 'datasets':
            show_datasets(args, known)
            return 0
        names = selected_datasets(args.datasets, known)
        if args.execute and args.command in ('pilot', 'full'):
            if not args.db.is_dir() or not args.virus_map.is_file():
                raise RuntimeError('数据库或病毒映射尚未建立；先按 ./metag.sh --help 完成共享参考准备。')
        failures = []
        for name in names:
            root = args.data_root / name
            print(f'\n[{name}] {args.command}: {root}', flush=True)
            try:
                if not root.is_dir():
                    raise RuntimeError(f'数据目录不存在：{root}')
                if args.command == 'prepare':
                    prepare(args, name, root, known[name])
                elif args.command == 'pilot':
                    pilot(args, root)
                else:
                    full(args, root)
                print(f'[{"完成" if args.execute else "计划"}] {name}', flush=True)
            except (RuntimeError, ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
                failures.append(name)
                print(f'[失败] {name}: {error}', file=sys.stderr, flush=True)
        print(f'\n已选择 {len(names)} 个数据集；失败 {len(failures)} 个。'
              + (f' 待处理：{", ".join(failures)}' if failures else ''), flush=True)
        return 1 if failures else 0
    except (RuntimeError, ValueError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
