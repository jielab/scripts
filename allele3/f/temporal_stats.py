"""Shared estimands for PC responses and measured allele-pair responses."""
from __future__ import annotations
import csv
import gzip
import hashlib
import math
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t

NA = {'', '.', 'NA', 'nan', 'None', 'unknown'}


def number(value):
    try: value = float(value)
    except (TypeError, ValueError): return None
    return value if math.isfinite(value) else None


def seed_for(seed, key):
    return (int(seed) + int(hashlib.sha256(str(key).encode()).hexdigest()[:8], 16)) % 2**32


def write_table(path, records, fields=None):
    records = list(records)
    fields = fields or list(dict.fromkeys(k for row in records for k in row)) or ['status', 'reason']
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in records:
            writer.writerow({k: 'NA' if row.get(k) is None else row[k] for k in fields})


def read_table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def age_bin(age, width, maximum):
    if age is None or age < 0 or age > maximum: return None
    return min(int(age // width), math.ceil(maximum / width) - 1)


def technology(row):
    raw = row.get('data_types', '').lower()
    if 'shotgun' in raw and ('capture' in raw or '1240' in raw): return 'mixed_shotgun_capture'
    if 'shotgun' in raw: return 'shotgun_label'
    if 'capture' in raw or '1240' in raw: return 'capture_label'
    return 'other_reported' if raw not in NA else 'unknown'


def bootstrap_mean(values, replicates, seed):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1: values = values[:, None]
    if len(values) < 2 or replicates < 2:
        return None, None, None, None
    rng = np.random.default_rng(seed)
    draws = np.empty((replicates, values.shape[1]))
    for start in range(0, replicates, 32):
        size = min(32, replicates - start)
        idx = rng.integers(len(values), size=(size, len(values)))
        draws[start:start + size] = values[idx].mean(axis=1)
    lo, hi = np.quantile(draws, [.025, .975], axis=0)
    return draws.std(axis=0, ddof=1), lo, hi, draws


def fixed_composition(values, times, strata, weights, epochs, minimum=20, bootstrap=200, seed=1):
    """Fixed weights over a declared set; missing strata yield NA, never reweight."""
    values = np.asarray(values, dtype=float)
    if values.ndim == 1: values = values[:, None]
    times = np.asarray(times); strata = np.asarray(strata)
    if not weights or any(w <= 0 for w in weights.values()) or not math.isclose(sum(weights.values()), 1., abs_tol=1e-10):
        raise ValueError('Fixed positive weights must sum to one')
    result = []
    for epoch in epochs:
        cells = {g: np.flatnonzero((times == epoch) & (strata == g) & np.isfinite(values).all(axis=1)) for g in weights}
        missing = [g for g, idx in cells.items() if len(idx) < minimum]
        record = dict(time_bin=int(epoch), n=sum(map(len, cells.values())), n_eff=None, estimate=None,
                      se=None, ci_low=None, ci_high=None, status='not_estimable' if missing else 'estimated',
                      reason='insufficient_common_support:' + ';'.join(missing) if missing else '')
        if not missing:
            record['n_eff'] = 1. / sum(weights[g]**2 / len(idx) for g, idx in cells.items())
            record['estimate'] = sum(weights[g] * values[idx].mean(axis=0) for g, idx in cells.items())
            draws = []
            for g, idx in cells.items():
                boot = bootstrap_mean(values[idx], bootstrap, seed_for(seed, (epoch, g)))[3]
                if boot is not None: draws.append(weights[g] * boot)
            if len(draws) == len(weights):
                aggregate = sum(draws)
                record['se'] = aggregate.std(axis=0, ddof=1)
                record['ci_low'], record['ci_high'] = np.quantile(aggregate, [.025, .975], axis=0)
        record['stratum_n'] = {g: len(idx) for g, idx in cells.items()}
        result.append(record)
    return result


def linear_model(samples, response, *, pcs=None, minimum=20, condition_limit=1e6):
    """Ancient-only HC3 association, complete cases; PC responses use pcs=None."""
    y = np.asarray(response, dtype=float)
    use_coverage = any(number(s.get('target_coverage')) is not None for s in samples)
    categorical = ['region', 'measurement_source']
    active = [i for i, s in enumerate(samples) if s.get('cohort') == 'ancient' and s.get('formal', True)
              and number(s.get('age_bp')) is not None and np.isfinite(y[i])
              and all(s.get(k, 'unknown') not in NA for k in categorical)
              and (not use_coverage or number(s.get('target_coverage')) is not None)
              and (pcs is None or s['sample_id'] in pcs)]
    base = dict(n=len(active), n_eff=len(active), estimate=None, se=None, ci_low=None, ci_high=None, p=None, q=None,
                status='not_estimable', reason='', model='linear_HC3_complete_case', units='response_per_1000_years_toward_present')
    if len(active) < minimum:
        return dict(base, reason='insufficient_complete_cases')
    ss = [samples[i] for i in active]
    names = ['intercept', 'toward_present_per_1000yr']
    columns = [np.ones(len(ss)), np.array([-s['age_bp']/1000 for s in ss])]
    for key in categorical:
        levels = sorted({s[key] for s in ss})
        for level in levels[1:]:
            columns.append(np.array([float(s[key] == level) for s in ss])); names.append(key + ':' + level)
    if use_coverage:
        columns.append(np.array([float(s['target_coverage']) for s in ss])); names.append('target_coverage')
    if pcs is not None:
        matrix = np.array([pcs[s['sample_id']] for s in ss])
        for k in range(matrix.shape[1]): columns.append(matrix[:, k]); names.append(f'PC{k+1}_AVG')
    # Constant nuisance covariates are recorded, not mistaken for information.
    dropped = [names[i] for i in range(2, len(names)) if np.std(columns[i]) < 1e-12]
    keep = [0, 1] + [i for i in range(2, len(names)) if names[i] not in dropped]
    X = np.column_stack([columns[i] for i in keep]); names = [names[i] for i in keep]
    for j in range(1, X.shape[1]):
        X[:, j] -= X[:, j].mean()
        if j > 1: X[:, j] /= X[:, j].std()
    rank = int(np.linalg.matrix_rank(X)); condition = float(np.linalg.cond(X))
    base.update(design_columns=';'.join(names), constant_covariates=';'.join(dropped), design_rank=rank,
                condition_number=condition if math.isfinite(condition) else None,
                included_sample_ids=';'.join(s['sample_id'] for s in ss), coverage_adjusted=use_coverage)
    if rank != X.shape[1] or condition > condition_limit or np.std(X[:, 1]) < 1e-12:
        return dict(base, reason='time_source_coverage_collinearity')
    residual_t = X[:, 1] - np.delete(X, 1, axis=1) @ np.linalg.lstsq(np.delete(X, 1, axis=1), X[:, 1], rcond=None)[0]
    vif = float(np.sum(X[:, 1]**2)/np.sum(residual_t**2)); base['time_vif'] = vif
    if vif > 1000: return dict(base, reason='unstable_time_support_vif_above_1000')
    if len(ss) <= X.shape[1] + 5: return dict(base, reason='insufficient_residual_degrees_of_freedom')
    projection = np.linalg.pinv(X); beta = projection @ y[active]
    leverage = np.einsum('ij,ji->i', X, projection)
    if np.max(leverage) >= .999: return dict(base, reason='high_leverage')
    residual = y[active] - X @ beta
    se = float(np.sqrt(np.sum((projection[1] * residual / (1-leverage))**2)))
    if se < 1e-14: return dict(base, estimate=float(beta[1]), reason='zero_residual_variance_no_inference')
    df = len(ss)-X.shape[1]; crit = float(student_t.ppf(.975, df))
    return dict(base, estimate=float(beta[1]), se=se, ci_low=float(beta[1]-crit*se), ci_high=float(beta[1]+crit*se),
                p=float(2*student_t.sf(abs(beta[1]/se), df)), status='estimated', reason='')


def adjust_fdr(records):
    indices = [i for i, r in enumerate(records) if number(r.get('p')) is not None]
    indices.sort(key=lambda i: records[i]['p'])
    last = 1.
    for rank, i in reversed(list(enumerate(indices, 1))):
        last = min(last, records[i]['p'] * len(indices) / rank)
        records[i]['q'] = last
    return records
