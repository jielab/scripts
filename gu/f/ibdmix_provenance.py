"""Compare IBDmix run contracts without treating resource relocation as new data.

pipeline_version in ibdmix.sh is the semantic contract for the workflow scripts:
changes to genotype preparation, calling or filtering MUST bump that version.
Script file timestamps and mask manifests remain audit information. Actual mask
BED checksums and all input stat records remain part of reuse validation.
"""
from __future__ import annotations
import argparse
import difflib
from pathlib import Path
import sys

# Recorded resource migration; do not collapse arbitrary paths or basenames.
RESOURCE_RELOCATIONS = {'/mnt/i/refGen/': '/mnt/e/refGen/'}
WORKFLOW_SCRIPTS = {'/mnt/d/scripts/gu/f/ibdmix.sh',
                    '/mnt/d/scripts/gu/f/ibdmix_workflow.py'}


def contract(text):
    # GNU stat's format in older run records emitted literal backslash-t.
    rows = [line.replace('\\t', '\t').split('\t') for line in text.splitlines() if line]
    refs = next((row[1].split() for row in rows if row[0] == 'refs'), [])
    versions = [row for row in rows if row[0] == 'pipeline_version']
    masks = [row for row in rows if row[0] == 'excluded_mask_sha256']
    units = {row[1] for row in masks if len(row) == 4}
    units.update(row[1] for row in rows if row[0] == 'mask_manifest_sha256' and len(row) >= 2)
    if len(versions) != 1 or len(versions[0]) != 2 or not versions[0][1] or not refs or not units:
        raise ValueError('Incomplete provenance: pipeline version, references and mask checksums are required')
    expected = {(unit, ref) for unit in units for ref in refs}
    actual = {(row[1], row[2]) for row in masks if len(row) == 4 and len(row[3]) == 64
              and all(c in '0123456789abcdef' for c in row[3])}
    if actual != expected or len(masks) != len(expected):
        raise ValueError('Incomplete or duplicate excluded-mask checksums')
    result = []
    for row in rows:
        key = row[0]
        if key in {'mask_root', 'mask_manifest_sha256'}:
            continue  # Actual BED hashes, rather than their storage metadata.
        if key == 'software' and len(row) == 2:
            path = row[1].rsplit(':', 2)[0]
            if path in WORKFLOW_SCRIPTS:
                continue  # Explicit pipeline_version controls semantic compatibility.
        for i, value in enumerate(row):
            for old, new in RESOURCE_RELOCATIONS.items():
                if value.startswith(old):
                    row[i] = new + value[len(old):]
        result.append('\t'.join(row))
    return sorted(result)


def differences(previous, current):
    before, after = contract(previous), contract(current)
    return ''.join(difflib.unified_diff([s+'\n' for s in before], [s+'\n' for s in after],
                                      fromfile='previous run contract', tofile='current run contract'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('previous', type=Path)
    parser.add_argument('current', type=Path)
    args = parser.parse_args()
    try:
        diff = differences(args.previous.read_text(), args.current.read_text())
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 1
    if diff:
        print(diff, end='', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
