#!/usr/bin/env python3
"""Prepare the audited GCST007020 source; this is not a general Y/Z decoder.

Y=first ID allele and Z=second ID allele was checked against 1000G EUR
frequencies in both orientations. Every recovered REF is checked against hg19.
The downloaded original is retained. Downstream use is through gwas_format.sh.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from gwas_alleles import IndexedFasta, id_alleles

SOURCE_SHA256 = '14d54ed2a090a86ea2e65b83d83489e4359e019c6e3b9cb5c2e0ab9526798190'
HEADER = 'SNP CHR POS EA NEA EAF N BETA SE P LOG10P'.split()


def decode(row, fasta):
    snp, c, pos, _, ea, nea, eaf, missing, beta, se, _, p = row
    reasons = []
    if c == '25':
        # PLINK 25 is XY; this study's source uses X coordinates for PAR.
        if not (60001 <= int(pos) <= 2699520 or 154931044 <= int(pos) <= 155260560):
            raise ValueError('Source CHR=25 outside GRCh37 X PAR: ' + snp)
        c = '23'
        reasons.append('PLINK_XY_25_to_X_23')
    if {ea, nea} == {'Y', 'Z'}:
        pair = id_alleles(snp, pos)
        if pair:
            if fasta.fetch(c, int(pos), len(pair[0])) != pair[0]:
                raise ValueError('ID reference allele disagrees with GRCh37: ' + snp)
            alleles = dict(zip(('Y', 'Z'), pair))
            ea, nea = alleles[ea], alleles[nea]
            reasons.append('Y_first_ID_allele_Z_second_ID_allele')
        else:
            ea = nea = 'NA'
            reasons.append('symbolic_SV_without_sequence_or_END;association_retained_alleles_NA')
    elif not set(ea + nea) <= set('ACGT'):
        raise ValueError('Unexpected allele coding: ' + snp)
    miss = float(missing)
    if not 0 <= miss < 1:
        raise ValueError('Invalid F_MISS: ' + snp)
    n = str(int(205327 * (1 - miss) + 0.5))
    lp = format(-math.log10(float(p)), '.12g') if float(p) > 0 else 'NA'
    return [snp, c, pos, ea, nea, eaf, n, beta, se, p, lp], reasons


def run(a):
    src, dst, qc = Path(a.input), Path(a.output), Path(a.qc_prefix)
    if src.resolve() == dst.resolve():
        raise ValueError('Preserve the original source')
    if dst.exists():
        raise ValueError('Prepared output already exists; inspect its audit before replacing it')
    dst.parent.mkdir(parents=True, exist_ok=True)
    qc.parent.mkdir(parents=True, exist_ok=True)
    fasta = IndexedFasta(a.fasta)
    counts, chromosomes = Counter(), Counter()
    digest = hashlib.sha256()
    with tempfile.TemporaryDirectory(prefix='mpb.prepare.', dir=dst.parent) as tmpdir:
        tsv = Path(tmpdir) / 'rows.tsv'
        audit = Path(tmpdir) / 'audit.tsv.gz'
        with src.open('rb') as inp, tsv.open('w') as out, gzip.open(audit, 'wt') as log:
            first = next(inp); digest.update(first)
            expected = 'SNP CHR BP GENPOS ALLELE1 ALLELE0 A1FREQ F_MISS BETA SE P_BOLT_LMM_INF P_BOLT_LMM'.split()
            if first.decode().split() != expected:
                raise ValueError('Unexpected source header')
            log.write('SNP\tCHR_SOURCE\tEA_SOURCE\tNEA_SOURCE\tCHR\tEA\tNEA\tREASON\n')
            for n, line in enumerate(inp, 1):
                digest.update(line)
                row = line.decode().rstrip('\r\n').split('\t')
                result, reasons = decode(row, fasta)
                out.write('\t'.join(result) + '\n')
                counts['input'] += 1
                chromosomes[result[1]] += 1
                for reason in reasons:
                    counts[reason] += 1
                if reasons:
                    log.write('\t'.join([row[0], row[1], row[4], row[5], result[1], result[3], result[4], ';'.join(reasons)]) + '\n')
                if n % 2000000 == 0:
                    print(f'Prepared {n:,} rows', flush=True)
        if digest.hexdigest() != SOURCE_SHA256:
            raise ValueError('Source differs from the audited GCST007020 download; no final output published')
        packed = Path(tmpdir) / 'prepared.gz'
        with packed.open('wb') as out:
            bg = subprocess.Popen(['bgzip', '-@', '4', '-c'], stdin=subprocess.PIPE, stdout=out)
            try:
                bg.stdin.write(('\t'.join(HEADER) + '\n').encode())
                with subprocess.Popen(['sort', '-T', tmpdir, '-S', '1G', '-t', '\t', '-k2,2n', '-k3,3n', '-k1,1', str(tsv)],
                                      stdout=subprocess.PIPE, env=dict(os.environ, LC_ALL='C')) as sorter:
                    shutil.copyfileobj(sorter.stdout, bg.stdin)
                    if sorter.wait():
                        raise RuntimeError('Sort failed')
                bg.stdin.close()
                if bg.wait():
                    raise RuntimeError('BGZF compression failed')
            except BaseException:
                bg.kill(); bg.wait()
                raise
        subprocess.run(['tabix', '-s', '2', '-b', '3', '-e', '3', '-S', '1', str(packed)], check=True)
        os.replace(packed, dst)
        os.replace(str(packed)+'.tbi', str(dst)+'.tbi')
        os.replace(audit, str(qc)+'.source.alleles.tsv.gz')
    fasta.close()
    Path(str(dst)+'.grch').write_text('37\n')
    report = dict(source=str(src), source_sha256=digest.hexdigest(), output=str(dst), counts=counts,
                  chromosomes=chromosomes, build=37, n_total=205327,
                  note='BETA, SE, EAF and P_BOLT_LMM preserved. Symbolic SV association rows retained with NA alleles.')
    Path(str(qc)+'.source.prepare.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ('input', 'output', 'qc-prefix'):
        p.add_argument('--'+arg, required=True)
    p.add_argument('--fasta', default='/mnt/e/refGen/fasta/GRCH37.fasta')
    run(p.parse_args())
