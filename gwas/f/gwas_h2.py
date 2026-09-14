#!/usr/bin/env python3
"""Summary-statistic SNP heritability, with explicit estimands and unavailable components."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from gwas_ldsc import DEFAULT_REF, DEFAULT_WEIGHTS, DEFAULT_ALLELES, validate_references, cache_paths
from gwas_liftover import identity


def parse_h2(text):
    pattern = r'Total Observed scale h2:\s*([-+\deE.]+)\s*\(([-+\deE.]+)\)'
    match = re.search(pattern, text)
    if not match:
        raise ValueError('No observed-scale h2 estimate found in LDSC log')
    import math
    estimate, se = map(float, match.groups())
    if not all(map(math.isfinite, (estimate, se))) or se < 0:
        raise ValueError('LDSC returned non-finite h2/SE')
    return estimate, se


def report_rows(trait, sex, estimate, se, chromosomes):
    auto = set(str(c) for c in range(1, 23))
    extra = set(chromosomes) - auto
    rows = []
    def add(metric, value='NA', error='NA', scope='all', status='UNAVAILABLE', reason=''):
        rows.append([trait, metric, value, error, 'LDSC' if value != 'NA' else 'not_estimated',
                     'observed', scope, sex, status, reason])
    add('h2_autosome', estimate, se, '1-22', 'ESTIMATED',
        'SNP heritability tagged by the specified LD reference; not pedigree heritability')
    if extra:
        add('h2_total', reason='Cannot call autosomal LDSC whole-genome h2; unestimated chromosomes: ' + ','.join(sorted(extra)))
    else:
        add('h2_total', estimate, se, '1-22', 'ESTIMATED', 'Input contains autosomes only; same estimate as h2_autosome')
    for name, code in [('X', '23'), ('Y', '24')]:
        add('h2_chr' + name, scope=name, status='UNAVAILABLE' if code in chromosomes else 'NOT_PRESENT',
            reason='Standard LDSC/reference does not support this chromosome; requires a validated sex/ploidy-specific LD method' if code in chromosomes else 'No input variants on this chromosome')
    add('h2_sig', reason='Do not run LDSC on GWAS-selected significant SNPs: selection bias. Requires an appropriate LD-aware variance/partition model and phenotype-scale information or individual data')
    add('h2_sig_chrX', scope='X', reason='Same limitation as h2_sig, plus male-X dosage/LD scaling; not the X proportion of significant-SNP heritability')
    add('h2_related', reason='Pedigree/close-relative variance requires individual phenotypes and family/genotype data; a reference kinship matrix alone is insufficient')
    for target in ('male', 'female'):
        if sex == target:
            add('h2_' + target + 's', estimate, se, '1-22', 'ESTIMATED',
                'Autosomal estimate for the declared sex-specific GWAS; not an independent estimate or all-chromosome GREML')
        else:
            add('h2_' + target + 's', reason='Requires a corresponding sex-specific GWAS or individual-level analysis; cannot split pooled summary statistics')
    return rows


def run(a):
    out = Path(a.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    done = out / 'h2.done.json'
    table = out / (a.trait + '.h2.tsv')
    refs = validate_references(a.ref_ld_chr, a.w_ld_chr)
    signature = dict(version=2, source=identity(a.gwas_file), sex=a.sex,
                     references=[identity(f) for f in refs], alleles=identity(a.merge_alleles),
                     python=a.python, conda_env=a.conda_env)
    if a.replace == 'FALSE' and done.exists() and table.exists():
        previous = json.loads(done.read_text())
        if previous.get('signature') == signature and previous.get('table') == identity(table):
            print('Reuse h2: ' + str(table)); return
    done.unlink(missing_ok=True)
    cmd = [sys.executable, str(Path(__file__).with_name('gwas_ldsc.py')),
           '--gwas-files', a.gwas_file, '--output-dir', str(out / 'ldsc'),
           '--merge-alleles', a.merge_alleles, '--ref-ld-chr', a.ref_ld_chr,
           '--w-ld-chr', a.w_ld_chr, '--run-rg', 'FALSE', '--conda-env', a.conda_env]
    if a.python:
        cmd += ['--python', a.python]
    with (out / 'h2.run.log').open('w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    stem = Path(a.gwas_file).name.removesuffix('.gz')
    estimate, se = parse_h2((out / 'ldsc' / 'h2.log' / (stem + '.h2.log')).read_text())
    with cache_paths(Path(a.gwas_file))[2].open() as handle:
        chromosomes={row['CHR']:int(row['INPUT']) for row in csv.DictReader(handle,delimiter='\t')}
    rows = report_rows(a.trait, a.sex, estimate, se, chromosomes)
    tmp = table.with_suffix('.tmp')
    with tmp.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['GWAS', 'METRIC', 'ESTIMATE', 'SE', 'METHOD', 'SCALE', 'CHROMOSOMES', 'SEX', 'STATUS', 'REASON'])
        writer.writerows(rows)
    os.replace(tmp, table)
    (out / 'h2.components.log').write_text('\n'.join('\t'.join(map(str, row)) for row in rows) + '\n')
    (out / 'h2.err').unlink(missing_ok=True)
    done.write_text(json.dumps(dict(signature=signature, table=identity(table),
                                   status='completed_supported_estimates; see component statuses'), indent=2))
    print(str(table))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gwas-file', required=True)
    p.add_argument('--trait', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--sex', choices=('unknown', 'mixed', 'male', 'female'), default='unknown')
    p.add_argument('--ref-ld-chr', default=DEFAULT_REF)
    p.add_argument('--w-ld-chr', default=DEFAULT_WEIGHTS)
    p.add_argument('--merge-alleles', default=DEFAULT_ALLELES)
    p.add_argument('--conda-env', default='ldsc')
    p.add_argument('--python')
    p.add_argument('--replace', choices=('TRUE', 'FALSE'), default='FALSE')
    args = p.parse_args()
    try:
        run(args)
    except Exception as exc:
        out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
        (out / 'h2.done.json').unlink(missing_ok=True)
        (out / 'h2.err').write_text(str(exc) + '\nSee h2.run.log for tool output.\n')
        raise
