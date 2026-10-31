"""Strict shared-coordinate contract; historical hashes are evidence, not guesses."""
from __future__ import annotations
import gzip
import json
import math
import os
from pathlib import Path
import re
from analysis_integrity import identity, sha256

DEFAULT_PROVENANCE = '/mnt/f/analysis/allele3_migration_20261028/all_files_verified.json'


def historical_identity(path, receipt=None, recorded=None):
    actual = identity(path, full=True)
    expected = (recorded or {}).get('sha256')
    provenance = None
    receipt = receipt or os.environ.get('ALLELE_PCA_PROVENANCE', DEFAULT_PROVENANCE)
    if not expected and receipt and Path(receipt).is_file():
        evidence = json.loads(Path(receipt).read_text())
        if evidence.get('status') != 'verified':
            raise ValueError('Historical provenance receipt is not verified')
        matches = [r for r in evidence.get('records', [])
                   if str(Path(r.get('destination', '')).resolve()) == actual['path']
                   or (recorded and str(Path(evidence.get('source', '')) / r.get('source_relative', '')) == recorded.get('path'))]
        if len(matches) > 1: raise ValueError('Ambiguous historical provenance entry')
        if matches:
            expected = matches[0].get('sha256'); provenance = identity(receipt, full=True)
    # The original metadata adapter saved full output hashes alongside the
    # metadata. This receipt survives removal of D outputs / old audit folders.
    adapter = Path(path).resolve().parent / 'pca3d_metadata.json'
    if not expected and adapter.is_file():
        saved_outputs = json.loads(adapter.read_text()).get('outputs', {})
        matches = [(name, digest) for name, digest in saved_outputs.items()
                   if name == actual['path'] or (recorded and name == recorded.get('path'))
                   or Path(name).name == Path(path).name]
        if len(matches) != 1:
            raise ValueError('Ambiguous/missing original metadata-adapter output receipt')
        expected = matches[0][1]
        if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-fA-F]{64}', expected):
            raise ValueError('Invalid original metadata-adapter digest')
        provenance = identity(adapter)
    if expected and actual['sha256'] != expected:
        raise ValueError(f'Historical content changed: {path}')
    return {'input': actual, 'historical_metadata_verification': 'verified' if expected else 'unverified',
            'expected_sha256': expected, 'provenance_receipt': provenance}


def read_scores(path, count, *, basis_dir=None, receipt=None, require_verified=True):
    if count < 0: raise ValueError('PC count must be nonnegative')
    if count == 0:
        return {}, [], {'expected_pc_count': 0, 'actual_pc_count': 0, 'coordinate_kind': 'disabled_by_pcs_0',
                        'basis_id': None, 'status': 'explicitly_unadjusted'}
    if not path: raise ValueError(f'{count} PCs requested but no coordinates supplied')
    path = Path(path).resolve()
    opener = gzip.open if path.suffix == '.gz' else open
    scores, records, row_bases = {}, {}, set()
    with opener(path, 'rt') as handle:
        lines = (line for line in handle if line.strip() and not line.startswith('##'))
        try: names = [x.lstrip('#') for x in next(lines).split()]
        except StopIteration: raise ValueError('Empty PCA coordinates')
        if len(set(names)) != len(names) or 'IID' not in names:
            raise ValueError('Coordinates need unique headers and IID')
        if 'PC1_AVG' not in names and any(re.fullmatch(r'PC\d+_SUM', x) for x in names):
            raise ValueError('PC SUM scores are not accepted; supply verified projected AVG on a common marker panel')
        if any(re.fullmatch(r'PC\d+', x) for x in names):
            raise ValueError('Raw eigenvectors are not projected AVG coordinates')
        columns = [f'PC{i}_AVG' for i in range(1, count + 1)]
        missing = set(columns) - set(names)
        if missing: raise ValueError('Missing consecutive requested PCs: ' + ','.join(sorted(missing)))
        for line in lines:
            cells = line.split()
            if len(cells) != len(names): raise ValueError('Malformed PCA coordinate row')
            row = dict(zip(names, cells)); iid = row['IID']
            if iid in {'', '.', 'NA'} or iid in scores: raise ValueError('Duplicate/missing PCA IID: ' + iid)
            values = [float(row[k]) for k in columns]
            if not all(math.isfinite(v) for v in values): raise ValueError('Nonfinite PC coordinate: ' + iid)
            if row.get('coordinate_kind', 'projected_average') not in {'projected_average', 'projected_AVG'}:
                raise ValueError('Mixed coordinate kinds')
            if row.get('basis_id'): row_bases.add(row['basis_id'])
            scores[iid] = values; records[iid] = row
    if not scores: raise ValueError('No coordinate rows')
    if len(row_bases) > 1: raise ValueError('Mixed PCA basis IDs')
    audit = {'expected_pc_count': count, 'actual_pc_count': len(columns), 'coordinate_kind': 'projected_average',
             'coordinates': identity(path), 'basis_id': None, 'status': 'unverified', 'n_coordinates': len(scores)}
    # Automatic discovery is restricted to a completed receipt in the coordinate's
    # own directory, never a guess based on column names or another default basis.
    # Producer receipts may only be auto-discovered beside their exact score file.
    if not receipt:
        candidate = path.with_name(path.stem + ".coordinate_receipt.json")
        if candidate.is_file():
            receipt = candidate
    basis = Path(basis_dir).resolve() if basis_dir else path.parent
    plan_path = basis / 'pca_plan.json'
    if plan_path.is_file() and path == basis / 'all.sscore':
        from pca_project import validate_completed_coordinates
        verified = validate_completed_coordinates(basis)
        plan = json.loads(plan_path.read_text())
        audit.update(status='verified', basis_id=verified['signature'], basis_dir=str(basis),
                     source_basis_id=verified['signature'],
                     coordinate_frame_id='original-frozen-coordinates:'+verified['signature'],
                     basis_receipts={p: identity(basis / p, full=True) for p in ('pca_plan.json', 'completed.sha256', 'PCA_COMPLETE')},
                     metadata_recorded=plan.get('inputs', {}).get('metadata'))
    elif receipt:
        proof = json.loads(Path(receipt).read_text())
        if proof.get('schema') == 'allele3-pca-coordinate-receipt-1':
            from pca_coordinate_receipts import verify_projection_receipt
            audit.update(verify_projection_receipt(path, receipt, columns, basis_dir))
        elif proof.get('schema'):
            raise ValueError('Unsupported coordinate receipt schema')
        else:
            # Keep the existing explicitly supplied legacy receipt contract.
            audit['coordinates'] = identity(path, full=True)
            if (proof.get('status') != 'verified' or proof.get('coordinate_kind') != 'projected_average'
                    or proof.get('coordinates_sha256') != audit['coordinates']['sha256'] or not proof.get('basis_id')
                    or not proof.get('marker_sha256')):
                raise ValueError('Coordinate receipt does not verify AVG, basis, marker set and full score SHA')
            audit.update(status='verified', basis_id=proof['basis_id'], supplied_receipt=identity(receipt, full=True),
                         metadata_recorded={'sha256': proof.get('metadata_sha256')} if proof.get('metadata_sha256') else None)
    elif basis_dir:
        raise ValueError('Explicit basis does not bind these coordinates; supply their projection receipt')
    if row_bases and (audit['basis_id'] not in row_bases): raise ValueError('Row basis disagrees with verified receipt')
    if require_verified and audit['status'] != 'verified':
        raise ValueError('No verified shared-basis coordinate receipt; adjusted model unavailable')
    return scores, columns, audit


def coordinate_join(samples, scores):
    ids = [s['sample_id'] for s in samples]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate metadata sample_id')
    return {'n_metadata': len(ids), 'n_coordinates': len(scores), 'n_joined': len(set(ids) & scores.keys()),
            'metadata_without_coordinates': sorted(set(ids) - scores.keys()),
            'coordinates_without_metadata': sorted(scores.keys() - set(ids))}
