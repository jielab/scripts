"""Reuse completed logical shards and bundle only unbuilt inputs, without gaps."""
import json
from pathlib import Path

from library import digest

BUNDLE_SHARDS = 4  # Four <=512 MiB inputs, while individual records stay ~1 MiB.
BUNDLE_VERSION = 'metacache-bundle-v1'


def bundle_identity(coverage):
    return digest([BUNDLE_VERSION, coverage, {'k': 16, 'sketch': 16, 'window': 127, 'stride': 112}])


def bundle_job(group, records, dbroot):
    first, last = records[0][0], records[-1][0]
    coverage = [{'group': group, 'shard': j, 'key': key} for j, _, key in records]
    return {'db': dbroot/f'part{group:04d}b{first:05d}-{last:05d}',
            'inputs': [folder for _, folder, _ in records],
            'identity': bundle_identity(coverage), 'coverage': coverage}


def completed_bundle(group, start, key, dbroot, complete):
    matches = []
    for count in range(1, BUNDLE_SHARDS + 1):
        records = [(j, None, key) for j in range(start, start + count)]
        job = bundle_job(group, records, dbroot)
        try:
            saved = json.loads(Path(str(job['db'])+'.complete.json').read_text())
            if saved.get('coverage') == job['coverage'] and complete(job['db'], job['identity']):
                matches.append((count, job))
        except (OSError, ValueError):
            pass
    if len(matches) > 1:
        raise RuntimeError('Overlapping completed bundles; refusing duplicate coverage')
    return matches[0] if matches else None


def plan_group(group, stream, dbroot, complete):
    """Consume each logical input exactly once. Existing small databases win.

    New receipts bind the physical DB to its precise ordered logical coverage.
    A final short bundle remains reusable even if preparation was interrupted.
    """
    pending = []
    records = iter(enumerate(stream, 1))
    for j, (folder, key) in records:
        legacy = dbroot/f'part{group:04d}s{j:05d}'
        bundled = completed_bundle(group, j, key, dbroot, complete)
        if complete(legacy, key):
            if bundled:
                raise RuntimeError('Completed small database overlaps a completed bundle')
            if pending:
                yield bundle_job(group, pending, dbroot)
                pending = []
            yield {'db': legacy, 'inputs': [folder], 'identity': key,
                   'coverage': [{'group': group, 'shard': j, 'key': key}]}
        elif bundled:
            if pending:
                yield bundle_job(group, pending, dbroot)
                pending = []
            count, job = bundled
            combined = [(j, folder, key)]
            for expected in range(j+1, j+count):
                item = next(records, None)
                if item is None or item[0] != expected or item[1][1] != key:
                    raise RuntimeError('Completed bundle coverage exceeds/differs from current inputs')
                if complete(dbroot/f'part{group:04d}s{expected:05d}', key):
                    raise RuntimeError('Completed bundle overlaps a completed small database')
                if completed_bundle(group, expected, key, dbroot, complete):
                    raise RuntimeError('Overlapping completed bundles')
                combined.append((expected, item[1][0], key))
            yield bundle_job(group, combined, dbroot)
        else:
            pending.append((j, folder, key))
            if len(pending) == BUNDLE_SHARDS:
                yield bundle_job(group, pending, dbroot)
                pending = []
    if pending:
        yield bundle_job(group, pending, dbroot)
