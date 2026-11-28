"""Derived summaries of frozen windows and chromosome-paired AVI profiles."""
import csv
from collections import Counter, defaultdict
from itertools import groupby

import numpy as np
from research_common import bh, chrom_key, number, open_text


def rows(path, required=()):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f'{path} lacks columns {sorted(missing)}')
        for row in reader:
            if None in row or None in row.values():
                raise ValueError(f'{path}: incomplete TSV row at line {reader.line_num}')
            yield row


def checked(row, name):
    value = number(row[name])
    if value is None:
        raise ValueError(f'Expected finite {name} in matched pair {row["chrom"]}:{row["pos"]}')
    return value


def control_set(row):
    ids = row['control_positions'].split(',')
    weights = list(map(float, row['control_weights'].split(',')))
    if (row['control_chrom'] != row['chrom'] or len(ids) != int(row['n_controls'])
            or len(set(ids)) != len(ids) or len(ids) != len(weights)
            or not all(np.isfinite(w) and w > 0 for w in weights)
            or not np.isclose(sum(weights), 1.0)):
        raise ValueError(f'Invalid matched control set at {row["chrom"]}:{row["pos"]}')
    return sorted(zip(ids, weights))


def near_far(path, near, far, bootstrap, seed, qualifying_alts=None):
    columns = ('chrom', 'pos', 'inner_bp', 'outer_bp', 'metric', 'case_value', 'control_mean',
               'control_chrom', 'control_positions', 'control_weights', 'n_controls')
    blocks = defaultdict(lambda: defaultdict(lambda: np.zeros(5)))
    counts = defaultdict(Counter)
    previous = None
    for (chrom, pos), records in groupby(rows(path, columns), lambda r: (r['chrom'], int(r['pos']))):
        site = chrom_key(chrom), pos
        if previous is not None and previous >= site:
            raise ValueError(f'Expected sorted, unique case groups in {path} at {chrom}:{pos}')
        previous = site
        selected = defaultdict(dict)
        for r in records:
            if qualifying_alts is not None and r.get('qualifying_alts') != str(qualifying_alts):
                continue
            band = (int(r['inner_bp']), int(r['outer_bp']))
            m = r['metric']
            selected[m]  # Retain cases with coverage only in intermediate distance bands.
            if band in (near, far):
                if band in selected[m]:
                    raise ValueError(f'Duplicate band for {m} at {chrom}:{pos}')
                selected[m][band] = r
        for m, bands in selected.items():
            count = counts[m]
            count['represented'] += 1
            if near not in bands or far not in bands:
                count['incomplete_bands'] += 1
                continue
            n, f = bands[near], bands[far]
            if control_set(n) != control_set(f):
                count['different_controls'] += 1
                continue
            blocks[m][chrom] += [checked(n,'case_value'), checked(f,'case_value'),
                                 checked(n,'control_mean'), checked(f,'control_mean'), 1]
    output = []
    for m in sorted(counts):
        by_chr = np.asarray(list(blocks[m].values()))
        result = dict(metric=m, near_inner_bp=near[0], near_outer_bp=near[1],
                      far_inner_bp=far[0], far_outer_bp=far[1], n_represented=counts[m]['represented'],
                      missing_a_band=counts[m]['incomplete_bands'], different_controls=counts[m]['different_controls'],
                      n_chromosomes=len(by_chr), n_cases=0, status='no_common_pairs', p=None, fdr=None,
                      case_near=None, case_far=None, control_near=None, control_far=None, case_change=None,
                      control_change=None, case_change_ci_low=None, case_change_ci_high=None,
                      control_change_ci_low=None, control_change_ci_high=None, estimate=None, ci_low=None, ci_high=None)
        if len(by_chr):
            total = by_chr.sum(axis=0)
            cn, cf, bn, bf = total[:4] / total[4]
            effect = (cn - cf) - (bn - bf)
            result.update(case_near=cn, case_far=cf, control_near=bn, control_far=bf, case_change=cn-cf,
                          control_change=bn-bf, estimate=effect, n_cases=int(total[4]), status='one_chromosome')
            if len(by_chr) >= 2:
                rng = np.random.default_rng(seed)
                draws = np.empty((bootstrap, 3))
                # Keep case weights intact. Each draw samples entire chromosomes.
                for begin in range(0, bootstrap, 256):
                    sample = rng.integers(len(by_chr), size=(min(256, bootstrap-begin), len(by_chr)))
                    s = by_chr[sample].sum(axis=1)
                    case_changes = (s[:,0]-s[:,1])/s[:,4]
                    control_changes = (s[:,2]-s[:,3])/s[:,4]
                    draws[begin:begin+len(s)] = np.column_stack((case_changes, control_changes, case_changes-control_changes))
                lo, hi = np.quantile(draws, [.025,.975], axis=0)
                p = (1+np.count_nonzero(abs(draws[:,2]-effect) >= abs(effect))) / (bootstrap+1)
                result.update(case_change_ci_low=lo[0], case_change_ci_high=hi[0], control_change_ci_low=lo[1],
                              control_change_ci_high=hi[1], ci_low=lo[2], ci_high=hi[2], p=p, status='chromosome_bootstrap')
        output.append(result)
    for row, q in zip(output, bh([r['p'] for r in output])):
        row['fdr'] = q
    return output


def window_summary(path, density):
    keys = ['threshold','window_bp'] + (['channel'] if density else [])
    groups = defaultdict(Counter)
    required = keys + (['mds_z'] if density else ['p_value','nominated'])
    for row in rows(path, required):
        group = groups[tuple(row[k] for k in keys)]
        group['n_windows'] += 1
        if density:
            z = number(row['mds_z'])
            if z is not None:
                group['n_scored'] += 1
                group['sum_z'] += z
        else:
            group['n_tested'] += number(row['p_value']) is not None
            group['n_nominated'] += int(row['nominated'])
    results = []
    for key, g in sorted(groups.items()):
        r = dict(zip(keys,key), n_windows=g['n_windows'])
        if density:
            r.update(n_scored=g['n_scored'], mean_scored_window_z=g['sum_z']/g['n_scored'] if g['n_scored'] else None)
        else:
            r.update(n_tested=g['n_tested'], n_nominated=g['n_nominated'])
        results.append(r)
    return results


def simulation_summary(path):
    groups = defaultdict(list)
    for row in rows(path, ['scenario','persistence_frequency','n_above_persistence_threshold','fixation_at_least_99pct']):
        key = (row['scenario'], row['persistence_frequency'])
        groups[key].append((int(row['n_above_persistence_threshold']),int(row['fixation_at_least_99pct'])))
    return [dict(scenario=s, persistence_threshold=f, replicates=len(v),
                 fraction_three_or_more=sum(n>=3 for n,_ in v)/len(v),
                 fraction_four=sum(n==4 for n,_ in v)/len(v),
                 fraction_fixed_99pct=sum(fix for _,fix in v)/len(v)) for (s,f),v in sorted(groups.items())]


def carrier_summary(path):
    keys = ('reference','family','group','allele_role','tract_call')
    groups = defaultdict(Counter)
    for row in rows(path, (*keys,'n_tested','n_called','n_carriers','AC','AN')):
        g = groups[tuple(row[k] for k in keys)]
        g['allele_rows'] += 1
        for k in ('n_tested','n_called','n_carriers','AC','AN'):
            g[k] += int(row[k])
    result = []
    for key,g in sorted(groups.items()):
        result.append(dict(zip(keys,key), n_allele_rows=g['allele_rows'], tested_person_allele_pairs=g['n_tested'],
                           called_person_allele_pairs=g['n_called'], carrier_person_allele_pairs=g['n_carriers'],
                           sum_AC=g['AC'],sum_AN=g['AN'], copy_weighted_frequency=g['AC']/g['AN'] if g['AN'] else None))
    return result


def time_models(path):
    # Some original fields list every individual in the model. Keep only
    # estimates and QC suitable for a manuscript table.
    fields = ('family','partition','role','membership_scope','region','measure','estimand','model',
              'units','n','n_eff','n_fixed_pairs','estimate','se','ci_low','ci_high','p','q',
              'status','reason','genetic_adjustment','inference_scope','n_chromosome_sensitivity',
              'chromosome_jackknife_se','chromosome_sensitivity_status')
    prior = csv.field_size_limit()
    try:
        csv.field_size_limit(max(prior,64*1024*1024))
        for row in rows(path,('estimate','q','membership_scope','role','n_fixed_pairs')):
            yield {key:row.get(key,'NA') for key in fields}
    finally:
        csv.field_size_limit(prior)


def fixed_pair_counts(path):
    fields = ('family','membership_scope','partition','role','status')
    counts = Counter(tuple(r[k] for k in fields) for r in rows(path,fields))
    for group,n in sorted(counts.items()):
        yield dict(zip(fields,group),n_pair_rows=n)
