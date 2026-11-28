"""Stage order and deterministic per-dataset batching for the public CLI."""
import argparse
from collections import Counter
import re

ORDER = ('group', 'species', 'evaluate', 'report')
PUBLIC = {'split-group': 'group', 'split-species': 'species', 'evaluate': 'evaluate', 'report': 'report'}


def actions(value):
    if value == 'all':
        return list(ORDER)
    if value == 'split':
        return ['group', 'species', 'report']
    result = [PUBLIC.get(x.strip(), x.strip()) for x in value.split(',')]
    if not result or any(x not in ORDER for x in result):
        raise argparse.ArgumentTypeError('Choose all or comma-separated split-group,split-species,evaluate,report')
    if len(set(result)) != len(result) or result != sorted(result, key=ORDER.index):
        raise argparse.ArgumentTypeError('Stages must be unique and ordered: split-group,split-species,evaluate,report')
    return result


def positive(value):
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError('Batch size must be a positive integer') from None
    if number < 1:
        raise argparse.ArgumentTypeError('Batch size must be a positive integer')
    return number


def select_batches(value, total, *, available_only=False):
    if value in (None, 'all'):
        return list(range(1, total + 1))
    selected = set()
    for item in value.split(','):
        item = item.strip()
        match = re.fullmatch(r'([1-9][0-9]*)(?:-([1-9][0-9]*)?)?', item)
        if not match:
            raise ValueError('Use --batch 1-2, 1,3-5, 3-, or all')
        first = int(match[1])
        last = (int(match[2]) if match[2] else (max(first, total) if available_only else total)) if '-' in item else first
        if first > last or (last > total and not available_only):
            raise ValueError(f'Batch selection {item!r} is outside 1-{total}')
        values = set(range(first, last + 1))
        if selected & values:
            raise ValueError('Batch selection contains repeated indices')
        selected.update(values)
    return sorted(n for n in selected if n <= total)


def cohort_batch_indices(samples, size, selection):
    """Select only available batches; a requested upper limit may exceed all datasets."""
    size = positive(size)
    counts = Counter(s['dataset'] for s in samples)
    total = max(((count + size - 1) // size for count in counts.values()), default=0)
    return select_batches(selection, total, available_only=True)


def batches(samples, size, selection):
    selected = cohort_batch_indices(samples, size, selection)
    result = []
    for dataset in sorted({s['dataset'] for s in samples}):
        members = sorted((s for s in samples if s['dataset'] == dataset), key=lambda s: s['sample_id'])
        chunks = [members[i:i + size] for i in range(0, len(members), size)]
        for number in selected:
            if number > len(chunks):
                continue
            result.append((dataset, number, chunks[number - 1]))
    return result
