"""Delete inventoried raw FASTQ only after durable, validated dataset completion."""
from datetime import datetime, timezone
import errno
import json
import os
from pathlib import Path
import stat
import tempfile
import time
import uuid


def plain_parents(path):
    for parent in path.parents:
        if parent.is_symlink():
            raise RuntimeError(f'raw 清理拒绝目录链接：{parent}')


def sync_directory(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        try:
            os.fsync(fd)
        except OSError as error:
            # Some mounted filesystems cannot fsync directories. File fsync is
            # still mandatory, including the completion marker and deletion log.
            if error.errno not in (errno.EINVAL, errno.ENOTSUP):
                raise
    finally:
        os.close(fd)


def durable_json(path, data):
    fd, tmp = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
        sync_directory(path.parent)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def source_identity(path, raw):
    plain_parents(path)
    info = path.lstat()
    target = str(path)
    target_info = info
    if stat.S_ISLNK(info.st_mode):
        target = str(path.resolve())
        target_info = path.stat()
    if not stat.S_ISREG(target_info.st_mode):
        raise RuntimeError(f'raw 清理对象不是普通 FASTQ：{path}')
    inventory = [str(path.relative_to(raw)), info.st_size, info.st_mtime_ns,
                 target, target_info.st_size, target_info.st_mtime_ns]
    return dict(inventory=inventory, inode=info.st_ino, device=info.st_dev,
                is_symlink=stat.S_ISLNK(info.st_mode))


def sync_outputs(root, state):
    """Flush already verified products without a second full FASTQ read."""
    files = {item['path']: item for item in state['artifacts']}
    for run in state['runs']:
        saved = json.loads((root / run['stem'] / 'split.done').read_text())
        files.update({item['path']: item for item in saved['outputs']})
    directories = {root}
    for name, expected in files.items():
        path = Path(name)
        if not path.is_relative_to(root) or '..' in path.parts or path.is_symlink():
            raise RuntimeError(f'raw 清理发现异常输出路径：{path}')
        plain_parents(path)
        with path.open('rb') as stream:
            info = os.fstat(stream.fileno())
            if (not stat.S_ISREG(info.st_mode) or info.st_size != expected['bytes']
                    or info.st_mtime_ns != expected['mtime_ns']):
                raise RuntimeError(f'raw 清理前输出发生变化：{path}')
            os.fsync(stream.fileno())
        directories.update(p for p in path.parents if p.is_relative_to(root))
    # Persist products before the dataset-level completion marker, then their
    # directory entries. Neither this marker nor any output is changed here.
    with (root / 'split.done').open('rb') as stream:
        os.fsync(stream.fileno())
    for folder in sorted(directories, key=lambda p: len(p.parts), reverse=True):
        sync_directory(folder)
    sync_directory(root.parent)


def remove_source(item, raw):
    path = Path(item['path'])
    if source_identity(path, raw) != item['identity']:
        raise RuntimeError(f'raw 文件在清理规划后发生变化，已停止：{path}')
    # unlink removes only this entry. A file symlink's target is never deleted.
    path.unlink()


def cleanup_dataset(root, raw, dataset, context):
    # Imported lazily because local_pipeline invokes this helper after splitting.
    from local_pipeline import NAME, dataset_receipt, digest, restore_dataset

    root, raw = Path(root).absolute(), Path(raw).absolute()
    source = raw / dataset
    if (not dataset or dataset in ('.', '..') or '/' in dataset or '\\' in dataset
            or root == raw or root.is_relative_to(raw) or raw.is_relative_to(root)):
        raise RuntimeError('raw 清理路径无效，输出必须与原始目录分离')
    plain_parents(root / 'split.done')
    if (root / 'split.done').is_symlink():
        raise RuntimeError('raw 清理拒绝链接形式的 split.done')
    state = dataset_receipt(root, raw, dataset)
    if state is None:
        raise RuntimeError(f'raw 清理要求数据集级 split.done：{root}')
    restore_dataset(root, raw, dataset, state, context)
    # Validate every candidate before deleting any. Removed originals are allowed
    # on resume; new or changed FASTQ are rejected by restore_dataset above.
    files, seen = [], set()
    for row in state['inventory']:
        relative = Path(row[0])
        if (relative.is_absolute() or '..' in relative.parts or not relative.parts
                or relative.parts[0] != dataset or not NAME.fullmatch(relative.name)
                or row[0] in seen):
            raise RuntimeError(f'raw 清单包含不安全或重复路径：{row[0]}')
        seen.add(row[0])
        path = raw / relative
        plain_parents(path)
        if not os.path.lexists(path):
            continue
        info = source_identity(path, raw)
        if info['inventory'] != list(row):
            raise RuntimeError(f'raw 文件与拆分完成清单不符：{path}')
        for suffix in ('.aria2', '.tmp'):
            if os.path.lexists(path.with_name(path.name + suffix)):
                raise RuntimeError(f'raw 文件仍有下载/临时标记：{path}')
        files.append(dict(path=str(path), identity=info))
    if not files:
        print(f'[metag] SKIP raw cleanup dataset={dataset}; 已处理原始 FASTQ 均已移除', flush=True)
        return dict(status='empty', deleted_files=0, deleted_bytes=0)
    print(f'[metag] START raw cleanup dataset={dataset}; files={len(files)}; '
          '验证输出并同步落盘后删除已处理 FASTQ', flush=True)
    sync_outputs(root, state)
    # Do not accept a completion receipt replaced while outputs were flushed.
    if dataset_receipt(root, raw, dataset) != state:
        raise RuntimeError('raw 清理前 split.done 发生变化')
    folder = root / 'raw_cleanup'
    if folder.is_symlink():
        raise RuntimeError(f'raw 清理记录目录不能是链接：{folder}')
    folder.mkdir(exist_ok=True)
    token = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '.' + uuid.uuid4().hex[:8]
    plan = folder / (token + '.plan.json')
    events = folder / (token + '.events.jsonl')
    report = folder / (token + '.summary.json')
    summary = dict(schema='metag-raw-cleanup-v1', dataset=dataset, raw_directory=str(source),
                   status='running', split_done_signature=digest(state), started_utc=datetime.now(timezone.utc).isoformat(),
                   planned_files=len(files), deleted_files=0, deleted_bytes=0, removed_symlinks=0,
                   plan=str(plan), events=str(events))
    durable_json(plan, dict(split_done_signature=digest(state), files=files))
    durable_json(report, summary)
    try:
        with events.open('x') as stream:
            last = time.monotonic()
            for item in files:
                remove_source(item, raw)
                summary['deleted_files'] += 1
                summary['deleted_bytes'] += item['identity']['inventory'][1]
                summary['removed_symlinks'] += int(item['identity']['is_symlink'])
                stream.write(json.dumps(dict(path=item['path'], status='deleted',
                                            symlink_target_preserved=item['identity']['is_symlink']),
                                        ensure_ascii=False) + '\n')
                stream.flush()
                os.fsync(stream.fileno())
                if time.monotonic() - last >= 60:
                    print(f'[metag] raw cleanup dataset={dataset}; {summary["deleted_files"]}/{len(files)} files', flush=True)
                    last = time.monotonic()
        for directory in {Path(item['path']).parent for item in files}:
            sync_directory(directory)
        summary['status'] = 'complete'
    except BaseException as error:
        summary.update(status='interrupted', error=str(error))
        raise
    finally:
        summary['finished_utc'] = datetime.now(timezone.utc).isoformat()
        durable_json(report, summary)
    durable_json(root / 'raw.cleanup.json', summary)
    print(f'[metag] DONE raw cleanup dataset={dataset}; deleted_files={summary["deleted_files"]}; '
          f'file_bytes={summary["deleted_bytes"]}; symlinks={summary["removed_symlinks"]}; receipt={report}', flush=True)
    return summary
