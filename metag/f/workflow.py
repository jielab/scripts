#!/usr/bin/env python3
"""One foreground workflow from downloaded inputs to per-dataset summaries."""
import argparse
import configparser
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from datasets import positive

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent / 'metag.sh'
VMR_NAME = 'VMR_MSL41.v1.20260729.xlsx'
VMR_SHA = 'b79b5d82a1b3b8e9dd5e19afe8fe1a8f441267474918a7cefa8ae4913adf45bb'


def raw_cleanup_default(path=None):
    config = configparser.ConfigParser(interpolation=None)
    path = path or HERE.parent / 'metag.conf'
    if Path(path).exists():
        with Path(path).open() as stream:
            config.read_file(stream)
    return config.getboolean('split', 'delete_raw_after_split', fallback=False)


def boolean_setting(value):
    if value.upper() not in ('TRUE', 'FALSE'):
        raise argparse.ArgumentTypeError('Use TRUE or FALSE')
    return value.upper() == 'TRUE'


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    try:
        delete_raw = raw_cleanup_default()
    except (OSError, ValueError, configparser.Error) as error:
        p.error(f'Invalid metag.conf: {error}')
    p.add_argument('--delete-raw-after-split', type=boolean_setting, default=delete_raw, metavar='TRUE|FALSE',
                   help='数据集拆分完成并校验后删除已处理原始 FASTQ [metag.conf，默认 FALSE]')
    p.add_argument('--dir-metag', type=Path, default=Path('/mnt/f/metag'))
    p.add_argument('--dir-annot', type=Path, default=Path('/mnt/f/annot/kraken2_ictv'))
    p.add_argument('--dir-fastq', type=Path, default=Path('/mnt/f/metag/raw'))
    p.add_argument('--dir-out', '--out', dest='out', type=Path, default=Path('/mnt/d/analysis/metag'))
    p.add_argument('--dir-split', type=Path, help='拆分 FASTQ 与逐数据集完成标记 [raw 同级的 split]')
    p.add_argument('--dat-list', default=None, help='默认 raw 下所有一级目录')
    p.add_argument('--stage', choices=['all', 'split'], default='all')
    p.add_argument('--classification-mode', choices=['quick', 'full'], default='full',
                   help='full = search both mates completely [default]; quick = early-stop screening')
    p.add_argument('--confidence', type=float, default=0.1)
    p.add_argument('--hit-groups', type=positive, default=3)
    p.add_argument('--memory-mapping', action='store_true')
    p.add_argument('--threads', type=positive, default=16, help='total Kraken2 thread budget; divided between split workers')
    # Accept old command lines, but these options no longer enable FASTQ QC.
    p.add_argument('--qc-jobs', type=positive, default=4, help=argparse.SUPPRESS)
    p.add_argument('--split-jobs', type=positive, default=4, help='parallel split runs [4]; enables shared memory mapping')
    p.add_argument('--qc-reads', type=int, default=10000, help=argparse.SUPPRESS)
    p.add_argument('--keep-assignments', action='store_true', help='write per-read assignment TSV.gz (off by default)')
    p.add_argument('--download-threads', type=positive, default=8)
    p.add_argument('--mask-jobs', type=positive, default=2,
                   help='parallel reference genomes, independent of download workers [2]')
    p.add_argument('--masker-threads', type=positive, default=1,
                   help='threads per reference masker [1]')
    p.add_argument('--prepare-memory-gb', type=positive, default=12,
                   help='reference download/masking RAM cap in GiB [12]; index build retains its own scope')
    p.add_argument('--max-db-gb', type=positive, default=24)
    p.add_argument('--dry-run', action='store_true', help='只显示顺序和命令，不写入或执行')
    p.add_argument('--reference-only', action='store_true', help='仅参考准备与审计（prep_ref 使用）')
    return p


def workflow(args):
    if not args.reference_only:
        from local_pipeline import run as local_run
        import copy
        reference_args = copy.copy(args)
        reference_args.reference_only = True
        return local_run(args, lambda: workflow(reference_args))
    ref = args.dir_annot.resolve()
    db = ref / '0.metadata/kraken2_broad_v1'
    map_dir = ref / '0.metadata/virus_map_v1'
    mapping = map_dir / 'virus_genome_types.sqlite'
    common = ['--db', db, '--virus-map', mapping, '--threads', args.threads]

    def run(*argv, capture=False):
        command = [str(x) for x in (ENTRY, *argv)]
        print(shlex.join(command), flush=True)
        if not args.dry_run:
            result = subprocess.run(command, check=True, text=True,
                                    stdout=subprocess.PIPE if capture else None)
            return result.stdout.strip() if capture else None

    # Reference preparation is required even when FASTQ contents are trusted.
    if not args.dry_run:
        if not ref.is_dir():
            raise RuntimeError(f'参考目录不存在：{ref}')
        markers = next(ref.rglob('*.aria2'), None)
        if markers:
            raise RuntimeError(f'参考仍在下载：{markers}；下载完成后重跑同一命令。')
    print('\n[metag] START reference preparation：恢复/复用参考库与病毒映射', flush=True)
    run('build-db', '--db', db, '--threads', args.threads,
        '--download-threads', args.download_threads, '--max-db-gb', args.max_db_gb,
        '--mask-jobs', min(args.mask_jobs, max(1, args.threads // args.masker_threads)),
        '--masker-threads', args.masker_threads, '--prepare-memory-gb', args.prepare_memory_gb, '--execute')
    if not mapping.is_file():
        if map_dir.exists() and not args.dry_run:
            raise RuntimeError(f'病毒映射目录不完整，保留现场，请检查：{map_dir}')
        vmr_dir = ref / '0.metadata/ictv_vmr_v1'
        url = f'https://ictv.global/sites/default/files/VMR/{VMR_NAME}'
        receipt = vmr_dir / VMR_NAME.replace('.xlsx', '.source.json')
        if receipt.is_file():
            url = json.loads(receipt.read_text())['source_url']
        vmr = run('virus', 'fetch', '--url', url, '--out-dir', vmr_dir, capture=True)
        if args.dry_run:
            vmr = vmr_dir / VMR_NAME
        elif hashlib.sha256(Path(vmr).read_bytes()).hexdigest() != VMR_SHA:
            raise RuntimeError('ICTV 文件与固定版本校验值不同，已停止。')
        cache = ref / '0.metadata/ictv_accessions_v1'
        run('virus', 'cache-accessions', '--vmr', vmr, '--cache-dir', cache)
        run('virus', 'build', '--db', db, '--vmr', vmr,
            '--accession-taxid', cache / 'accession_taxid.tsv', '--out', map_dir)
    print('\n[metag] START reference/audit：正式索引及病毒 DNA/RNA 覆盖审计', flush=True)
    run('db-audit', *common)
    run('reference-clean', '--db', db, '--virus-map', mapping, '--execute')
    print('\n' + ('计划已显示，未执行。' if args.dry_run else
                  '[metag] DONE reference preparation and audit'), flush=True)


def main():
    args = parser().parse_args()
    try:
        if args.masker_threads > args.threads:
            raise ValueError('--masker-threads must not exceed the total --threads budget')
        if args.classification_mode == 'quick' and args.hit_groups < 2:
            raise ValueError('Kraken2 2.17.1 quick 模式要求 --hit-groups >= 2；单命中请用 --classification-mode full')
        if not 0 <= args.confidence <= 1:
            raise ValueError('--confidence 必须在 0 到 1 之间')
        if args.dry_run:
            workflow(args)
            return 0
        # One scope covers reference preparation and splitting, including children.
        if os.environ.get('METAG_WORKFLOW_SCOPED') != '1':
            subprocess.run(['systemctl', '--user', 'show', '--property=Version', '--value'],
                           check=True, stdout=subprocess.DEVNULL)
            env = dict(os.environ, METAG_WORKFLOW_SCOPED='1')
            return subprocess.run([
                'systemd-run', '--user', '--scope', '--quiet', '--collect',
                f'--unit=metag-workflow-{os.getpid()}-{time.time_ns()}',
                '--property=MemoryHigh=30G', '--property=MemoryMax=36G',
                '--property=MemorySwapMax=4G', '--property=OOMPolicy=stop',
                '--property=KillMode=control-group', sys.executable, str(Path(__file__)),
                *sys.argv[1:]], env=env).returncode
        ref = args.dir_annot.resolve()
        if not ref.is_dir():
            raise RuntimeError(f'参考目录不存在：{ref}')
        (ref / '0.metadata').mkdir(exist_ok=True)
        # Share the installer preparation lock; build-db additionally owns its build lock.
        with (ref / '0.metadata/.metag_prepare.lock').open('a') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise RuntimeError('同一参考目录已有准备/分析流程在运行，请等待完成。')
            workflow(args)
        return 0
    except (OSError, RuntimeError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f'流程停止：{error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
