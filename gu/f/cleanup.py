#!/usr/bin/env python3
"""Remove inactive GU scratch caches; results and source genotypes are never targets."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil

CACHE_ROOTS = tuple(Path('/tmp') / name for name in ('gu-ukb-prepared', 'gu-intermediate', 'gu-native-view'))


def process_references():
    # The finishing wrapper and its ancestors are expected to mention their own
    # work root. Other invocations, workers, readers and open files protect it.
    ancestors = set()
    pid = os.getpid()
    while pid > 1:
        ancestors.add(pid)
        try:
            status = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
            pid = int(status[1])
        except (FileNotFoundError, ProcessLookupError):
            break
    references = []
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit() or int(proc.name) in ancestors:
            continue
        try:
            # Only the current user's jobs can own these caches. Unknown access
            # to another user's processes does not imply that their files are ours.
            if proc.stat().st_uid != os.getuid():
                continue
            if (proc / 'stat').read_text().rsplit(')', 1)[1].split()[0] == 'Z':
                continue  # exited processes cannot retain cache descriptors
            if (proc / 'comm').read_text().strip() == '(sd-pam)':
                continue  # systemd's credential helper has protected /proc data
            for name in ('cmdline', 'environ', 'maps'):
                references.append((proc / name).read_bytes().decode(errors='replace'))
            for link in [proc / 'cwd', *(proc / 'fd').iterdir()]:
                try:
                    references.append(os.readlink(link))
                except (FileNotFoundError, ProcessLookupError):
                    pass
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError as exc:
            raise RuntimeError(f'Cannot inspect own process {proc.name}; refusing cleanup') from exc
    return references


def key(path):
    return hashlib.sha256(str(path).encode()).hexdigest()[:16]


def candidates(analysis):
    """Map generated scratch directories back to their durable result scope."""
    intermediate, views = {}, {}
    for work in analysis.iterdir():
        if not work.is_dir() or work.is_symlink():
            continue
        views[key(work)] = work
        for method in work.iterdir():
            if method.is_dir() and not method.is_symlink():
                for run in method.iterdir():
                    if run.is_dir() and not run.is_symlink():
                        intermediate[key(run)] = (work, run)
    for cache_root in CACHE_ROOTS:
        if not cache_root.is_dir() or cache_root.is_symlink():
            continue
        for path in cache_root.iterdir():
            if not path.is_dir() or path.is_symlink() or path.stat().st_uid != os.getuid():
                continue
            work = None
            if cache_root.name == 'gu-native-view':
                work = views.get(path.name)
            elif cache_root.name == 'gu-intermediate':
                mapped = intermediate.get(path.name)
                if mapped:
                    work, run = mapped
                    # Masks are disposable only once native results are archived.
                    if not any(run.glob('*.raw.tar.gz')):
                        continue
            else:
                roots = set()
                for env in path.glob('*/*/*/gu-target.env'):
                    for line in env.read_text().splitlines():
                        if line.startswith('export GU_ANALYSIS_ROOT='):
                            roots.add(Path(shlex.split(line.split('=', 1)[1])[0]).resolve())
                if len(roots) == 1:
                    work = roots.pop()
            if work is not None and work.parent == analysis and work.is_dir():
                yield path, work


def cleanup(analysis, work=None, apply=False):
    analysis = Path(analysis).resolve()
    work = Path(work).resolve() if work else None
    refs = process_references()
    rows = []
    for path, owner in candidates(analysis):
        if work is not None and owner != work:
            continue
        if any(str(path) in ref or str(owner) in ref for ref in refs):
            rows.append(dict(path=str(path), work=str(owner), status='active-retained'))
            continue
        size = sum(p.lstat().st_blocks * 512 for p in path.rglob('*') if not p.is_symlink())
        row = dict(path=str(path), work=str(owner), allocated_bytes=size, status='candidate')
        if apply:
            # Refresh immediately before deletion: a new reader may have started
            # since planning. Never unlink a busy cache on the basis of age alone.
            latest = process_references()
            if any(str(path) in ref or str(owner) in ref for ref in latest):
                row['status'] = 'active-retained'
            else:
                shutil.rmtree(path)
                row['status'] = 'deleted'
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root', default='/mnt/f/analysis/gu')
    parser.add_argument('--work', help='Only clean this completed result scope')
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--report')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()
    rows = cleanup(args.analysis_root, args.work, args.apply)
    report = dict(entries=rows, deleted_bytes=sum(r.get('allocated_bytes', 0) for r in rows if r['status']=='deleted'))
    if args.report:
        Path(args.report).write_text(json.dumps(report, indent=2)+'\n')
    if args.quiet:
        print(f"GU CLEANUP: released {report['deleted_bytes'] / 2**30:.2f} GiB; active scopes retained")
    else:
        print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
