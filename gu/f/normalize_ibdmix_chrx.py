#!/usr/bin/env python3
"""Normalize published male-X filenames and repair native read links.

Dry-run by default. --apply requires --receipt and refuses name collisions.
Internal caller IDs and data contents are preserved.
"""
import argparse
import json
import os
from pathlib import Path
from ibdmix_layout import chrx_output_name


def identity(path):
    s = path.stat()
    return [s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns]


def inventory(root):
    moves = []
    links = []
    for directory, dirs, names in os.walk(root):
        d = Path(directory)
        # Native autosomal runs cannot reference male-X artifacts.
        if d.name == 'ibdmix' and d.parent.name == 'tmp':
            dirs[:] = [n for n in dirs if n.startswith('chrX')]
        for name in names:
            p = d / name
            is_x = 'X_MALE' in name or 'chrX' in name or any(n.startswith('chrX') for n in p.relative_to(root).parts[:-1])
            if not is_x:
                continue
            if p.is_symlink():
                links.append((p, os.readlink(p)))
            elif 'X_MALE' in name:
                # Only published files; native IDs remain available to callers.
                if 'ibdmix' in p.relative_to(root).parts[1:-1]:
                    raise ValueError(f'Unpublished native file: {p}')
                target = p.with_name(chrx_output_name(name))
                if target != p:
                    moves.append(dict(old=str(p), new=str(target), identity=identity(p)))
    return moves, links


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    moves, links = inventory(args.root.resolve())
    mapping = {row['old']: row['new'] for row in moves}
    targets = set(mapping.values())
    assert len(targets) == len(moves), 'Conflicting destinations'
    for target in targets:
        assert not os.path.lexists(target), f'Destination exists: {target}'
    rewrites = []
    broken = []
    for p, old in links:
        target = Path(os.path.abspath(p.parent / old))
        newtarget = Path(mapping.get(str(target), str(target)))
        # Some older logs lost the chrX prefix while their read links did not.
        if 'X_MALE' in target.name and not target.exists():
            candidate = target.with_name(chrx_output_name(target.name))
            if str(candidate) in targets or candidate.exists():
                newtarget = candidate
        if not p.exists():
            broken.append(str(p))
        if newtarget != target:
            assert str(newtarget) in targets or newtarget.exists(), newtarget
            rewrites.append(dict(path=str(p), before=old, after=os.path.relpath(newtarget, p.parent)))
    plan = dict(status='planned', moves=moves, links=rewrites, broken_before=broken)
    print(f'Preflight: {len(moves)} renames, {len(rewrites)} link updates, {len(broken)} pre-existing broken X links', flush=True)
    if args.receipt:
        with args.receipt.open('x') as f:
            json.dump(plan, f, indent=2)
    if not args.apply:
        return
    assert args.receipt, '--apply requires --receipt'
    for row in moves:
        old, new = Path(row['old']), Path(row['new'])
        assert identity(old) == row['identity'], old
        assert not os.path.lexists(new), new
        os.rename(old, new)
    for row in rewrites:
        p = Path(row['path'])
        assert os.readlink(p) == row['before'], p
        temporary = p.with_name(p.name + '.chrx-link')
        assert not os.path.lexists(temporary), temporary
        temporary.symlink_to(row['after'])
        os.replace(temporary, p)
    for row in moves:
        assert identity(Path(row['new'])) == row['identity'], row
        assert not os.path.lexists(row['old']), row
    broken_after = [str(p) for p, _ in links if not p.exists()]
    assert not (set(broken_after) - set(broken)), broken_after
    for row in rewrites:
        assert Path(row['path']).exists(), row
    plan.update(status='verified', broken_after=broken_after)
    args.receipt.write_text(json.dumps(plan, indent=2) + '\n')
    print(f'Verified: {len(moves)} files retain inode, size and mtime; {len(rewrites)} links updated; {len(broken_after)} broken X links remain', flush=True)


if __name__ == '__main__':
    main()
