#!/usr/bin/env python3
"""Validate REGENIE sample-sex metadata and chrX labels without reading genotypes."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path


def phenotype_keys(path):
    with open(path) as stream:
        header = stream.readline().split()
        fi, ii = header.index('FID'), header.index('IID')
        return {(row[fi], row[ii]) for line in stream if (row := line.split())}


def sample_sexes(path, fam=False):
    result = {}
    with open(path) as stream:
        if fam:
            rows = stream
            fi, ii, si = 0, 1, 4
        else:
            header = stream.readline().split()
            header = [x.lstrip('#') for x in header]
            if 'IID' not in header or 'SEX' not in header:
                raise ValueError('chrX PSAM must contain IID and SEX (PLINK: 1=male, 2=female)')
            fi = header.index('FID') if 'FID' in header else None
            ii, si = header.index('IID'), header.index('SEX')
            rows = stream
        for line in rows:
            row = line.split()
            if not row:
                continue
            key = (row[fi] if fi is not None else '0', row[ii])
            if key in result:
                raise ValueError('Duplicate sample keys in genotype metadata')
            result[key] = row[si]
    return result


def audit(phe, fam, psam=None, pvar=None, build='37', sex='all', keep=None):
    selected = phenotype_keys(phe)
    if keep:
        allowed = set()
        with open(keep) as stream:
            for line in stream:
                row = line.split()
                if not row or row[0].lstrip('#') in ('FID', 'IID'):
                    continue
                if len(row) != 2:
                    raise ValueError('REGENIE --keep audit requires two columns: FID IID')
                allowed.add(tuple(row))
        selected &= allowed
    typed = sample_sexes(fam, fam=True)
    typed_selected = selected & typed.keys()
    if not typed_selected:
        raise ValueError('No phenotype samples match the Step 1 FAM')
    xsex = sample_sexes(psam) if psam else None
    candidates = typed_selected & xsex.keys() if xsex is not None else typed_selected
    if not candidates:
        raise ValueError('No complete-case samples overlap Step 1 and chrX')
    if psam or sex != 'all':
        unknown = sum(typed[k] not in ('1', '2') or (xsex is not None and xsex[k] not in ('1', '2')) for k in candidates)
        if unknown:
            raise ValueError(f'{unknown} candidate samples have unknown/non-PLINK sex; fix FAM/PSAM metadata explicitly')
        mismatches = sum(typed[k] != xsex[k] for k in candidates) if xsex is not None else 0
        if mismatches:
            raise ValueError(f'{mismatches} candidate samples have discordant sex in Step 1 FAM and chrX PSAM')
    if sex != 'all':
        code = '1' if sex == 'male' else '2'
        candidates = {k for k in candidates if typed[k] == code}
    if not candidates:
        raise ValueError('No analysis samples remain after the requested genotype-sex filter')
    result = {'build': build, 'sex_subset': sex, 'phenotype_complete_cases': len(selected),
              'step1_overlap': len(typed_selected), 'audited_overlap': len(candidates),
              'sex_counts': dict(Counter(typed[k] for k in candidates)),
              'genotype_encoding_checked': False,
              'dosage_model': 'REGENIE default: male non-PAR X on 0-2 scale; no recoding'}
    if pvar:
        low, high, last = {'37': (2699520, 154931044, 155270560),
                           '38': (2781479, 155701383, 156040895)}[build]
        labels, regions = Counter(), Counter()
        with open(pvar) as stream:
            for line in stream:
                if line.startswith('##'):
                    continue
                if line.startswith('#'):
                    header = line.split()
                    ci, pi = header.index('#CHROM'), header.index('POS')
                    break
            else:
                raise ValueError('Missing PVAR #CHROM/POS header')
            for line in stream:
                if not line.strip():
                    continue
                row = line.split(); chrom = row[ci].upper(); pos = int(row[pi])
                if chrom not in ('X', '23', 'XY', '25', 'PAR1', 'PAR2'):
                    raise ValueError('chrX PVAR contains non-X/PAR chromosome labels; do not combine Y with X')
                if not 0 < pos <= last:
                    raise ValueError('chrX position outside the selected genome-build chromosome length')
                region = 'PAR1' if pos <= low else 'PAR2' if pos >= high else 'nonPAR'
                if chrom in ('XY', '25', 'PAR1', 'PAR2') and (region == 'nonPAR' or (chrom in ('PAR1','PAR2') and chrom != region)):
                    raise ValueError('PVAR PAR label disagrees with the selected build; check --grch and genotype preparation')
                labels[chrom] += 1; regions[region] += 1
        if not labels:
            raise ValueError('No chrX variants found')
        result.update(chromosome_labels=dict(labels), regions=dict(regions), par_option='b'+build)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phe', required=True); p.add_argument('--fam', required=True)
    p.add_argument('--psam'); p.add_argument('--pvar'); p.add_argument('--keep')
    p.add_argument('--build', choices=['37','38'], required=True)
    p.add_argument('--sex', choices=['all','male','female'], default='all')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    if bool(a.psam) != bool(a.pvar): p.error('--psam and --pvar must be supplied together')
    report = audit(a.phe, a.fam, a.psam, a.pvar, a.build, a.sex, a.keep)
    Path(a.output).write_text(json.dumps(report, indent=2)+'\n')
    print(f'Genotype metadata audit passed: {report["audited_overlap"]} samples; sex={a.sex}; build={a.build}')

if __name__ == '__main__':
    main()
