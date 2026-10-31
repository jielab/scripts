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
