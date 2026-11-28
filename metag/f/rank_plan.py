"""Bound resident-index batches by both sample count and compressed input size."""
from collections import defaultdict
from pathlib import Path


def automatic_batches(samples, cfg):
    maximum = max(1, int(cfg.get('rank_auto_batch_samples', 10)))
    byte_limit = max(1, int(float(cfg.get('rank_auto_batch_gib', 100)) * (1 << 30)))
    grouped = defaultdict(list)
    for sample in samples:
        grouped[sample['dataset']].append(sample)
    plans = []
    for dataset, members in sorted(grouped.items()):
        batch, size, index = [], 0, 1
        for sample in sorted(members, key=lambda x: x['sample_id']):
            total = sum(Path(sample[key]).stat().st_size for key in ('r1', 'r2') if sample.get(key))
            if batch and (len(batch) >= maximum or size + total > byte_limit):
                plans.append((dataset, index, batch))
                batch, size, index = [], 0, index + 1
            batch.append(sample)
            size += total
        if batch:
            plans.append((dataset, index, batch))
    # Interleave datasets so a large dataset does not monopolize all three lanes.
    return sorted(plans, key=lambda item: (item[1], item[0]))
