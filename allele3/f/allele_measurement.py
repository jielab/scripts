#!/usr/bin/env python3
"""Shared nucleotide-count reader for allele3 FINAL.

PVAR row shape alone does not prove
that its reference dosage is an exact nucleotide count after multiallelic split.
No inversion of one split ALT is allowed. Native multiallelic PGEN is supported.
An encoded ALT is measured on the panel's hardcalled sample denominator; it is
not asserted to be a conditional ALT-versus-reference genotype. Biallelic/split
REF counts remain unavailable because these inputs carry no native-export
provenance proving that every non-target allele was preserved.
"""
from __future__ import annotations

# Allele Measurement
# ======================================================================
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Mapping, Sequence
import numpy as np

BASES = frozenset('ACGT')
MISSING = -9

@dataclass(frozen=True)
class PanelRecord:
    chrom: str
    pos: int
    index: int
    alleles: tuple[str, ...]  # actual PVAR order: REF, ALT1, ALT2, ...
    pass_record: bool
    hq_alts: frozenset[str]
    source_id: str = '.'

    def __post_init__(self):
        if self.pos < 1 or self.index < 0:
            raise ValueError('Invalid PVAR coordinate/index')
        if not 2 <= len(self.alleles) <= 4 or len(set(self.alleles)) != len(self.alleles):
            raise ValueError('Expected distinct SNV nucleotide labels')
        if not set(self.alleles) <= BASES:
            raise ValueError('Only canonical SNVs are accepted')
        if not self.hq_alts <= set(self.alleles[1:]):
            raise ValueError('HQ ALT identity is not in its PVAR row')


def dosage_plan(records: Sequence[PanelRecord], counted_base: str) -> dict:
    """Plan one actual nucleotide; unavailable remains unavailable, never zero.

    Direct ALT reads need no completeness assertion: they describe the encoded
    allele. REF reads require a native multiallelic PGEN row. A set of biallelic
    rows is not a provenance proof: even a singleton may be an incomplete split.
    """
    unavailable = dict(status='not_represented_in_panel', counted_base=counted_base)
    if counted_base not in BASES:
        raise ValueError('Counted nucleotide must be A/C/G/T')
    records = list(records)
    if not records:
        return unavailable
    if len({(r.chrom, r.pos) for r in records}) != 1:
        raise ValueError('One genomic locus is required')
    if len({r.index for r in records}) != len(records):
        raise ValueError('Duplicate PVAR row index')
    if len({r.alleles[0] for r in records}) != 1:
        return dict(unavailable, status='inconsistent_panel_ref')
    native = [r for r in records if len(r.alleles) > 2]
    if native:
        if len(records) != 1:
            return dict(unavailable, status='mixed_native_and_split_representation')
        r = native[0]
        if counted_base not in r.alleles:
            return unavailable
        if not r.pass_record:
            return dict(unavailable, status='pvar_filter_failed')
        i = r.alleles.index(counted_base)
        hq = counted_base in r.hq_alts if i else set(r.alleles[1:]) <= r.hq_alts
        if not hq:
            return dict(unavailable, status='not_in_high_quality')
        return dict(status='assayed_pass', mode='pgen_allele', pvar_index=r.index,
                    allele_idx=i, counted_base=counted_base, panel_alleles=list(r.alleles),
                    source_id=r.source_id,
                    measurement_scope='native_joint_nucleotide_count',
                    reference_count_scope='direct_native_joint_GT', invert_dosage=0)
    # Every row below is a biallelic representation. Unique ALT counts are usable
    # independently, but that says nothing about unlisted bases/REF complements.
    alt_rows = [r for r in records if r.alleles[1] == counted_base]
    if len(alt_rows) > 1:
        return dict(unavailable, status='duplicate_counted_alt_records')
    if alt_rows:
        r = alt_rows[0]
        if not r.pass_record:
            return dict(unavailable, status='pvar_filter_failed')
        if counted_base not in r.hq_alts:
            return dict(unavailable, status='not_in_high_quality')
        return dict(status='assayed_pass', mode='pgen_allele', pvar_index=r.index,
                    allele_idx=1, counted_base=counted_base, panel_alleles=list(r.alleles),
                    source_id=r.source_id,
                    measurement_scope='panel_encoded_ALT_count',
                    reference_count_scope='not_inferred', invert_dosage=0)
    ref = records[0].alleles[0]
    if counted_base != ref:
        return unavailable
    return dict(unavailable, status='reference_requires_native_joint_genotype')


def read_counted_dosage(reader, plan: Mapping, sample_count: int, audit: Counter | None = None):
    """Return diploid counts 0/1/2 or -9. All consumers must call this function."""
    if plan.get('status') != 'assayed_pass':
        raise ValueError('Cannot read an unavailable measurement plan')
    if plan.get('invert_dosage', 0):
        raise ValueError('Single-ALT dosage inversion is prohibited')
    if sample_count < 1:
        raise ValueError('sample_count must be positive')
    if reader.get_raw_sample_ct() != sample_count:
        raise ValueError('PGEN sample count differs from PSAM')
    stats = audit if audit is not None else Counter()
    n_variants = reader.get_variant_ct()
    def one(index, allele_idx):
        if not isinstance(index, int) or not 0 <= index < n_variants:
            raise ValueError('PVAR row index lies outside PGEN')
        if not isinstance(allele_idx, int) or not 0 <= allele_idx <= 3:
            raise ValueError('Invalid nucleotide allele index')
        value = np.empty(sample_count, dtype=np.int8)
        reader.read(index, value, allele_idx=allele_idx)
        if np.any(~np.isin(value, [MISSING, 0, 1, 2])):
            raise ValueError('Unexpected non-diploid PGEN hardcall')
        return value
    if plan['mode'] == 'pgen_allele':
        result = one(plan['pvar_index'], plan['allele_idx'])
    else:
        raise ValueError('Unknown nucleotide-count mode')
    stats['called_measurements'] += int((result >= 0).sum())
    stats['missing_measurements'] += int((result < 0).sum())
    return result


def panel_pvar(prefix):
    """Require one unambiguous companion PVAR, including its native allele list."""
    present = [p for p in (Path(str(prefix)+'.pvar'), Path(str(prefix)+'.pvar.zst')) if p.is_file()]
    if len(present) != 1:
        raise ValueError('Exactly one PVAR/PVAR.ZST matching the PGEN is required')
    return present[0]


@contextmanager
def open_pgen_with_pvar(pgen_path, pvar_path, sample_count):
    """Supply true allele offsets, required for native multiallelic PGEN rows."""
    import pgenlib
    pvar = pgenlib.PvarReader(os.fsencode(pvar_path))
    try:
        reader = pgenlib.PgenReader(os.fsencode(pgen_path), raw_sample_ct=sample_count, pvar=pvar)
        try:
            yield reader
        finally:
            reader.close()
    finally:
        pvar.close()


def iter_pgen_chunks(path, sample_count, selected, chunk_size, *, pvar_path=None, audit=None):
    """Stream planned nucleotide counts with bounded genotype memory."""
    if chunk_size < 1:
        raise ValueError('chunk_size must be positive')
    if pvar_path is None:
        path = Path(path)
        if path.suffix != '.pgen':
            raise ValueError('Pass pvar_path explicitly for a nonstandard PGEN filename')
        pvar_path = panel_pvar(str(path)[:-5])
    with open_pgen_with_pvar(path, pvar_path, sample_count) as reader:
        for start in range(0, len(selected), chunk_size):
            chunk = selected[start:start+chunk_size]
            values = np.empty((len(chunk), sample_count), dtype=np.int8)
            for index, entry in enumerate(chunk):
                values[index] = read_counted_dosage(reader, entry['measurement_plan'], sample_count, audit)
            yield start, values, reader.get_variant_ct()


def iter_primary_counts(prefix, sample_count, primary, *, audit=None):
    """Stream one prespecified primary nucleotide per modern target locus."""
    with open_pgen_with_pvar(str(prefix)+'.pgen', panel_pvar(prefix), sample_count) as reader:
        for entry in primary:
            yield entry, read_counted_dosage(reader, entry['measurement_plan'], sample_count, audit)


def select_primary_measurements(mapped, selected_alts, allele_audit):
    """Preserve a prespecified primary nucleotide; never substitute another ALT."""
    by_site = {}
    for entry in selected_alts:
        by_site.setdefault(entry['site_id'], {})[entry['allele']] = entry
    audit_lookup = {(r['site_id'], r['allele']): r for r in allele_audit}
    chosen, audit = [], []
    for target in mapped:
        site = target['site_id']
        wanted = target.get('primary_alt')
        explicit = wanted not in {None, '', '.', 'NA', 'nan', 'None'}
        if explicit and wanted not in target['alts'].split(','):
            raise ValueError(f'{site}: primary_alt must be one of the supplied modern ALTs')
        if not explicit:
            wanted = next(iter(sorted(by_site.get(site, {}))), None)
        entry = by_site.get(site, {}).get(wanted)
        status = audit_lookup.get((site, wanted), {}).get('panel_status', 'no_available_counted_base')
        rule = 'predefined_modern_primary_nucleotide' if explicit else 'deterministic_available_fallback_no_primary'
        row = dict(site_id=site, predefined_primary_alt38=target.get('primary_alt'),
                   primary_rule=rule, counted_base38=wanted,
                   primary_status='available' if entry is not None else status,
                   measurement_status=status)
        if entry is not None:
            chosen.append(dict(entry, primary_rule=rule, counted_base38=entry['allele']))
            row.update({k: entry[k] for k in ('counted_base37', 'pvar_index', 'allele_idx',
                                             'measurement_mode', 'counting_scope', 'invert_dosage')})
        audit.append(row)
    return chosen, audit


def phase_relation(records, sample):
    """VCF implicit phase sets: all phased calls without PS share one set."""
    if not all(bool(getattr(r.samples[sample], 'phased', False)) for r in records):
        return False, 'unphased_or_mixed'
    ps = [r.samples[sample].get('PS') for r in records]
    if all(x is None for x in ps):
        return True, 'implicit_phase_set'
    if all(x is not None for x in ps) and len(set(ps)) == 1:
        return True, 'explicit_phase_set'
    return False, 'different_or_mixed_phase_sets'


def index_target_records(targets, pvar_path, hq_path, *, rows, canonical_chrom):
    """Preserve full PVAR allele lists and every same-position companion row.

    ``rows`` is herv_time.rows; do not restrict to modern target ALTs before
    recording the representation. A hidden companion is precisely the REF risk.
    """
    wanted = {(t['chrom37'], int(t['pos37'])) for t in targets if t.get('mapping_status') == 'mapped'}
    hq = set()
    for row in rows(hq_path, ('CHROM','POS','REF','ALT'), whitespace=True):
        coord = canonical_chrom(row['CHROM']), int(row['POS'])
        if coord not in wanted or row.get('FILTER','PASS') not in {'PASS','.'}:
            continue
        ref = row['REF'].upper()
        for alt in row['ALT'].upper().split(','):
            if ref in BASES and alt in BASES and alt != ref:
                hq.add((*coord, ref, alt))
    found = {coord: [] for coord in wanted}
    invalid = {}
    for index,row in enumerate(rows(pvar_path, ('CHROM','POS','REF','ALT'), whitespace=True)):
        coord = canonical_chrom(row['CHROM']), int(row['POS'])
        if coord not in wanted:
            continue
        labels = (row['REF'].upper(), *row['ALT'].upper().split(','))
        # Keep a locus-level unavailable reason, never drop only the companion.
        if len(set(labels)) != len(labels) or not set(labels) <= BASES:
            invalid[coord] = 'non_snv_or_invalid_companion_at_locus'
            continue
        high = frozenset(a for a in labels[1:] if (*coord,labels[0],a) in hq)
        found[coord].append(PanelRecord(*coord,index,labels,row.get('FILTER','PASS') in {'PASS','.'},high,row.get('ID','.')))
    return {'records_by_position': found, 'invalid_by_position': invalid}


def catalog_measurements(mapped, indexed, *, orient, alt_role):
    """Return ALT plans, complete audit, and separate true-REF plans.

    REF curves are separate and are never duplicate ALT hypothesis tests.
    A changed build REF does not discard directly measurable other ALT bases.
    """
    selected, audit, reference = [], [], []
    common_fields = ('site_id','chrom','pos','ref','family','group','partition','comparison_partition',
                     'match_id','strict_filter_complete','mapping_status','locus_n_alt_ac5',
                     'locus_n_alt_ac10','locus_ac10_supported')
    for target in mapped:
        for base in [*target['alts'].split(','),target['ref']]:
            is_ref = base == target['ref']
            ac,an = target.get('_counts',{}).get(base,(None,None))
            entry = {k:target.get(k) for k in common_fields}
            entry.update(allele=base,role='REF' if is_ref else alt_role(target,base),modern_ac=ac,
                         modern_an=an,modern_af=ac/an if an else None,
                         modern_ac5=int(ac is not None and ac>=5),modern_ac10=int(ac is not None and ac>=10),
                         allele_id=f"{target['site_id']}:{base}",
                         contrast_id=f"{target['site_id']}|count:{base}",
                         comparison_kind='one_nucleotide_vs_all_other_nucleotides',
                         panel_status='mapping_unavailable',reference_count_scope='unavailable',invert_dosage=0)
            if target.get('mapping_status') == 'mapped':
                coordinate = target['chrom37'],int(target['pos37'])
                records = indexed['records_by_position'].get(coordinate,[])
                invalid_reason = indexed['invalid_by_position'].get(coordinate)
                base37 = orient(base,target['strand37to38'])
                entry['counted_base37'] = base37
                if invalid_reason:
                    plan = dict(status=invalid_reason)
                elif any(r.alleles[0] != target['ref37'] for r in records):
                    plan = dict(status='panel_reference_mismatch')
                else:
                    plan = dosage_plan(records,base37)
                entry['panel_status'] = plan['status']
                if plan['status'] == 'assayed_pass':
                    entry.update(measurement_plan=plan,pvar_index=plan['pvar_index'],allele_idx=plan['allele_idx'],
                                 measurement_mode=plan['mode'],counting_scope=plan['measurement_scope'],
                                 reference_count_scope=plan['reference_count_scope'],panel_ref37=target['ref37'],
                                 panel_alt37=','.join(plan['panel_alleles'][1:]),
                                 panel_alleles37=','.join(plan['panel_alleles']),source_id=plan['source_id'],
                                 variant_id=plan['source_id'] if plan['source_id'] != '.' else
                                     f"{coordinate[0]}_{coordinate[1]}_{','.join(plan['panel_alleles'])}")
                    if is_ref:
                        entry['reference_count_scope'] = ('direct_native_joint_modern_REF'
                            if plan['measurement_scope'] == 'native_joint_nucleotide_count'
                            else 'direct_encoded_ALT_is_modern_REF')
                    (reference if is_ref else selected).append(entry)
            audit.append(entry)
    for collection in (selected,reference):
        collection.sort(key=lambda r:(r['pvar_index'],r['allele']))
    if len({r['allele_id'] for r in selected+reference}) != len(selected)+len(reference):
        raise ValueError('Duplicate nucleotide measurement identity')
    return selected,audit,reference


# Temporal Stats
# ======================================================================
_description_temporal_stats = 'Shared estimands for PC responses and measured allele-pair responses.'

import csv
import gzip
import hashlib
import math

STATS_NA = {'', '.', 'NA', 'nan', 'None', 'unknown'}


def stats_number(value):
    try: value = float(value)
    except (TypeError, ValueError): return None
    return value if math.isfinite(value) else None


def stats_seed_for(seed, key):
    return (int(seed) + int(hashlib.sha256(str(key).encode()).hexdigest()[:8], 16)) % 2**32


def stats_write_table(path, records, fields=None):
    records = list(records)
    fields = fields or list(dict.fromkeys(k for row in records for k in row)) or ['status', 'reason']
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in records:
            writer.writerow({k: 'NA' if row.get(k) is None else row[k] for k in fields})


def stats_read_table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def stats_age_bin(age, width, maximum):
    if age is None or age < 0 or age > maximum: return None
    return min(int(age // width), math.ceil(maximum / width) - 1)


def stats_technology(row):
    raw = row.get('data_types', '').lower()
    if 'shotgun' in raw and ('capture' in raw or '1240' in raw): return 'mixed_shotgun_capture'
    if 'shotgun' in raw: return 'shotgun_label'
    if 'capture' in raw or '1240' in raw: return 'capture_label'
    return 'other_reported' if raw not in STATS_NA else 'unknown'


def stats_bootstrap_mean(values, replicates, seed):
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


def stats_fixed_composition(values, times, strata, weights, epochs, minimum=20, bootstrap=200, seed=1):
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
                boot = stats_bootstrap_mean(values[idx], bootstrap, stats_seed_for(seed, (epoch, g)))[3]
                if boot is not None: draws.append(weights[g] * boot)
            if len(draws) == len(weights):
                aggregate = sum(draws)
                record['se'] = aggregate.std(axis=0, ddof=1)
                record['ci_low'], record['ci_high'] = np.quantile(aggregate, [.025, .975], axis=0)
        record['stratum_n'] = {g: len(idx) for g, idx in cells.items()}
        result.append(record)
    return result


def stats_linear_model(samples, response, *, pcs=None, minimum=20, condition_limit=1e6):
    """Ancient-only HC3 association, complete cases; PC responses use pcs=None."""
    from scipy.stats import t as student_t
    y = np.asarray(response, dtype=float)
    use_coverage = any(stats_number(s.get('target_coverage')) is not None for s in samples)
    categorical = ['region', 'measurement_source']
    active = [i for i, s in enumerate(samples) if s.get('cohort') == 'ancient' and s.get('formal', True)
              and stats_number(s.get('age_bp')) is not None and np.isfinite(y[i])
              and all(s.get(k, 'unknown') not in STATS_NA for k in categorical)
              and (not use_coverage or stats_number(s.get('target_coverage')) is not None)
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


def stats_adjust_fdr(records):
    indices = [i for i, r in enumerate(records) if stats_number(r.get('p')) is not None]
    indices.sort(key=lambda i: records[i]['p'])
    last = 1.
    for rank, i in reversed(list(enumerate(indices, 1))):
        last = min(last, records[i]['p'] * len(indices) / rank)
        records[i]['q'] = last
    return records



