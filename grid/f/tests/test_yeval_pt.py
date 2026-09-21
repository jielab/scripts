"""PLINK integration: COJO allele/position matching, duplicate rsIDs and score sums."""
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


def make_genotypes(root, chromosome, records):
    vcf = root / f'input{chromosome}.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.2\n'
        f'##contig=<ID={chromosome}>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t1\t2\t3\n'
        + ''.join(f'{chromosome}\t{pos}\t{sid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t'
                  + '\t'.join(gt) + '\n' for pos, sid, ref, alt, gt in records))
    subprocess.run(['plink2', '--vcf', str(vcf), '--make-pgen', '--out',
                    str(root / f'gen/chr{chromosome}')], check=True, stdout=subprocess.DEVNULL)


with tempfile.TemporaryDirectory(prefix='yeval-pt-test-') as scratch:
    root = Path(scratch)
    (root / 'gen').mkdir()
    make_genotypes(root, 1, [
        (100, 'rs1', 'A', 'G', ['0/0', '0/1', '1/1']),
        (300, 'rs_at_pos', 'A', 'G', ['0/0', '0/0', '0/0']),
        (400, 'rs_at_pos', 'A', 'G', ['1/1', '0/1', '0/0']),
        (500, 'rs_wrong_pos', 'A', 'G', ['0/0', '0/1', '1/1']),
        (600, 'rs_bad_allele', 'A', 'C', ['0/0', '0/1', '1/1']),
        # A shared REF allele cannot distinguish these two split records.
        (700, 'rs_ambig', 'C', 'G', ['0/0', '0/1', '1/1']),
        (700, 'rs_ambig', 'C', 'T', ['0/0', '1/1', '0/1']),
        (800, 'rs_identical', 'A', 'G', ['0/0', '0/1', '1/1']),
        (800, 'rs_identical', 'A', 'G', ['0/0', '0/1', '1/1']),
        # LDL COJO contains full sequence alleles, not just single-base SNPs.
        (900, 'rs_insertion', 'C', 'CAA', ['0/0', '0/1', '1/1']),
        (1000, 'rs_deletion', 'TC', 'T', ['0/0', '0/1', '1/1']),
        (1100, 'rs_indel_mismatch', 'C', 'CA', ['0/0', '0/1', '1/1']),
    ])
    # Reproduce the real height failure, with different dosages for each ALT.
    make_genotypes(root, 3, [
        (13831614, 'rs12629061', 'C', 'G', ['1/1', '0/1', '0/0']),
        (13831614, 'rs12629061', 'C', 'T', ['0/0', '1/1', '0/1']),
    ])
    original_genotypes = {p: p.read_bytes() for p in (root / 'gen').glob('*.p*')}
    for j, pop in enumerate(['EUR', 'AFR', 'EAS', 'SAS']):
        source = root / 'gwas' / f'test.{pop}' / 'gwas'
        source.mkdir(parents=True)
        rows = [
            (1, 'rs1', 100, 'A' if j % 2 else 'G', j + 1),
            (1, 'rs_missing', 200, 'G', 100),
            (1, 'rs_at_pos', 400, 'G', 5),
            (1, 'rs_wrong_pos', 501, 'G', 100),
            (1, 'rs_bad_allele', 600, 'G', 100),
            (1, 'rs_ambig', 700, 'C', 100),
            (1, 'rs_identical', 800, 'G', 100),
            (1, 'rs_insertion', 900, 'CAA', 2),
            (1, 'rs_deletion', 1000, 'TC', 3),
            (1, 'rs_indel_mismatch', 1100, 'CAA', 100),
        ]
        if pop != 'SAS':
            rows.append((3, 'rs12629061', 13831614, ['T', 'G', 'C'][j], 3))
        pd.DataFrame(rows, columns=['Chr', 'SNP', 'bp', 'refA', 'bJ']).to_csv(
            source / f'test.{pop}.jma.cojo', sep='\t', index=False)
    cmd = [sys.executable, str(Path(__file__).resolve().parents[1] / 'yeval_pt.py'),
           '--trait', 'test', '--dir-gwas', str(root / 'gwas'),
           '--dir-gen', str(root / 'gen'), '--output', str(root / 'pt.pgs.gz'),
           '--remove', '']
    scored = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert 'PLINK command:' in scored.stdout and 'PLINK v' in scored.stdout
    scores = pd.read_csv(root / 'pt.pgs.gz', sep='\t')
    expected = dict(EUR=[16, 17, 9], AFR=[26, 15, 4], EAS=[16, 13, 10], SAS=[24, 14, 4])
    for pop, values in expected.items():
        np.testing.assert_allclose(scores[f'pt.{pop}'], values)
    matches = pd.read_csv(root / 'pt.pgs.gz.matches.tsv', sep='\t')
    eur = matches[matches.population == 'EUR'].set_index('SNP')
    assert eur.loc['rs12629061', 'selected_variant'] == '3:13831614:C:T'
    assert eur.loc['rs_at_pos', 'selected_variant'] == '1:400:A:G'
    assert eur.loc['rs_insertion', 'selected_variant'] == '1:900:C:CAA'
    assert eur.loc['rs_deletion', 'selected_variant'] == '1:1000:TC:T'
    for sid, status in [('rs_missing', 'absent'), ('rs_wrong_pos', 'position_mismatch'),
                        ('rs_bad_allele', 'allele_mismatch'), ('rs_ambig', 'ambiguous'),
                        ('rs_identical', 'ambiguous'), ('rs_indel_mismatch', 'allele_mismatch')]:
        assert eur.loc[sid, 'status'] == status
    audit = pd.read_csv(root / 'pt.pgs.gz.variants.tsv', sep='\t')
    assert (audit.requested == audit[['scored', 'absent', 'position_mismatch',
                                     'allele_mismatch', 'ambiguous']].sum(axis=1)).all()
    assert audit.duplicates_resolved.sum() == 6
    assert audit.ambiguous.sum() == 9
    for path, original in original_genotypes.items():
        assert path.read_bytes() == original, f'Source genotype was changed: {path}'
    stamp = (root / 'pt.pgs.gz').stat().st_mtime_ns
    sidecars = {p:p.read_bytes() for p in root.glob('pt.pgs.gz.*')}
    reused = subprocess.run(cmd + ['--out', str(root / 'cached-evaluation')],
                            check=True, capture_output=True, text=True)
    assert 'reuse matching cache' in reused.stdout
    assert (root / 'pt.pgs.gz').stat().st_mtime_ns == stamp
    assert not (root / 'cached-evaluation').exists()
    assert not list(root.glob('pt.*.log')) and not list(root.glob('*.jsonl'))
    for path, original in sidecars.items():
        assert path.read_bytes() == original
    # Invalid symbols/effects must still fail before touching the score cache.
    source = root / 'gwas/test.EUR/gwas/test.EUR.jma.cojo'
    weights = pd.read_csv(source, sep='\t')
    for column, value, message in [('refA','N','Invalid refA'), ('refA','<DEL>','Invalid refA'),
                                    ('refA','A,C','Invalid refA'), ('bJ',np.inf,'Invalid/nonfinite bJ'),
                                    ('bJ','invalid','Invalid/nonfinite bJ')]:
        bad = weights.astype({column:object})
        bad.loc[0,column] = value
        bad.to_csv(source,sep='\t',index=False)
        failed = subprocess.run(cmd,capture_output=True,text=True)
        assert failed.returncode != 0 and message in failed.stderr and 'rs1' in failed.stderr
        assert (root / 'pt.pgs.gz').stat().st_mtime_ns == stamp
    print('PT SNP/indel matching, duplicate rsIDs, score sums, invalid-input rejection and cache without duplicate outputs PASS')
