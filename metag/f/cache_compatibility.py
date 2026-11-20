"""Narrow, audited migration of classifications across non-classifier changes."""
from pathlib import Path
import json


def compatible_request(previous, current):
    if previous == current:
        return True
    if not isinstance(previous, dict) or not isinstance(current, dict):
        return False
    old, new = dict(previous), dict(current)
    old_code, new_code = old.pop('code', None), new.pop('code', None)
    if old != new or not old_code or not new_code:
        return False
    try:
        receipt = json.loads(Path(__file__).with_name('classification_compatibility.json').read_text())
    except (OSError, ValueError):
        return False
    return (receipt.get('current_code') == new_code and old_code in receipt.get('compatible_codes', []))


def legacy_agreement_valid(qc, summary, context, result, current_code):
    """Verify an old completed agreement against its exact sampled read digest."""
    import csv
    import random
    from collections import Counter, defaultdict
    from fastq_stream import fragments, new_digest, update_digest
    if (result.get('status') != 'complete' or result.get('full_database') != context['database']
            or result.get('matched_parameters') != context['parameters']
            or not compatible_request(summary.get('request'), dict(summary.get('request', {}), code=current_code))):
        return False
    qc = Path(qc)
    try:
        for name, stamp in context['inputs'].items():
            # Original classification observations include every sampled input.
            suffix = str(qc.relative_to(qc.parents[1])) + '/' + name
            expected = summary.get('artifact_observations', {}).get(suffix)
            if not expected or any(stamp.get(k) != v for k,v in expected.items()):
                return False
        with (qc/'sample.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        strata = defaultdict(list)
        for i,row in enumerate(rows):
            strata[row['stratum']].append(i)
        rng = random.Random(int(context['seed']) + 1)
        keep = set()
        counts = Counter()
        for key, indices in sorted(strata.items()):
            chosen = rng.sample(indices, min(int(context['pairs']), len(indices)))
            keep.update(chosen); counts[key] = len(chosen)
        selected = [row for i,row in enumerate(rows) if i in keep]
        digest = new_digest()
        for i,(a,b) in enumerate(fragments(qc/'sample.R1.fastq.gz',qc/'sample.R2.fastq.gz')):
            if i in keep:
                update_digest(digest,a,b)
        if result.get('query', {}).get('input_digest') != digest.hexdigest() or result.get('pairs') != len(selected):
            return False
        with (qc/'agreement/pairs.tsv').open(newline='') as handle:
            detail = list(csv.DictReader(handle, delimiter='\t'))
        if len(detail) != len(selected):
            return False
        for row, saved in zip(selected, detail):
            if (int(row['ordinal']) != int(saved['ordinal']) or row['group'] != saved['stage_group']
                    or int(row['taxid']) != int(saved['stage_taxid'])
                    or abs(float(saved['weight']) - int(row['stratum_pairs'])/counts[row['stratum']]) > 1e-9):
                return False
        return (qc/'agreement/transitions.tsv').is_file()
    except (OSError, ValueError, KeyError):
        return False
