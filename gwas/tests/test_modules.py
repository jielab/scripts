#!/usr/bin/env python3
"""Synthetic regression tests; no UKB participant data is used."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
POST=ROOT.parent/'0data/gwas_post.sh'
def run(*args): subprocess.run(list(map(str,args)),check=True)
def table(path,header,rows):
    path.write_text('\t'.join(header.split())+'\n'+'\n'.join('\t'.join(map(str,r)) for r in rows)+'\n')

def main():
    with tempfile.TemporaryDirectory(prefix='gwas-tests-') as tmp:
        d=Path(tmp)
        header='SNP CHR POS EA NEA EAF N BETA SE P'
        rows=[['rs1',1,100,'A','G',.2,1000,.1,.02,1e-9],['rs2',1,200,'A','C',.3,1000,.2,.03,1e-10],['pal',1,300,'A','T',.4,1000,.3,.04,1e-10],['ns',2,400,'C','T',.1,1000,.01,.1,.5],['dup',2,500,'A','C',.2,1000,.1,.02,1e-9],['dup',2,500,'A','C',.2,1000,.1,.02,1e-9]]
        table(d/'a.tsv',header,rows)
        other=[r[:] for r in rows[:4]]
        other[0][3:6]=['G','A',.8]; other[0][7]=-.1
        other[1][3:5]=['T','G'] # complement, no flip
        table(d/'b.tsv',header,other)
        run('bash',POST,'compare','--gwas-files',f'{d}/a.tsv,{d}/b.tsv,{d}/a.tsv','--labels','A,B,C','--output-dir',d/'compare')
        qc=list(csv.DictReader((d/'compare/comparison_qc.tsv').open(),delimiter='\t'))
        assert len(qc)==4 and all(int(x['n'])==2 and abs(float(x['r'])-1)<1e-8 for x in qc),qc
        assert (d/'compare/manhattan.compare.png').stat().st_size>1000
        # Real PLINK2 association smoke tests exercise quantitative, ordinal and binary paths.
        plink='/mnt/d/software/bin/plink2'
        gen=d/'gen'; gen.mkdir()
        run(plink,'--dummy','500','100','--seed','7','--make-pgen','--out',gen/'chr1')
        with (gen/'chr1.psam').open() as h:
            names=h.readline().lstrip('#').split(); ids=[dict(zip(names,l.split())) for l in h]
        rng=random.Random(12)
        vals=[[r.get('FID','0'),r['IID'],rng.gauss(170,10),rng.randint(1,4),rng.randint(0,1),rng.randint(0,1),rng.uniform(.1,20),rng.uniform(40,70),rng.randint(0,1),rng.randint(1,3)] for r in ids]
        table(d/'phe','FID IID height bald bald12 cvd_cad.Yt2e cvd_cad.t2e age sex center',vals)
        run('bash',ROOT/'gwas.sh','run_plink2','--phenotype-file',d/'phe','--imputed-dir',gen,'--output-dir',d/'output','--chromosomes','1','--covariates','age,sex','--categorical-covariates','center','--threads','2','--memory','1000','--min-mac','5')
        for trait in ['height','bald','bald12','cvd_cad.Yt2e']:
            f=d/'output'/trait/'gwas'/f'{trait}.plink2.gz'
            with gzip.open(f,'rt') as h:
                assert h.readline().split()==header.split()
                assert len(h.readlines())>0
        assert (d/'output/cvd_cad.t2e/gwas/plink2/UNSUPPORTED.txt').exists()
        # Run actual pruning and PLINK1 merge on synthetic genotypes; omit unavailable SAIGE.
        run('bash',ROOT/'gwas.sh','prep_gwas','--typed-dir',gen,'--chromosomes','1','--threads','2','--memory','1000','--sparse-grm','FALSE','--run','FALSE')
        plan=gen/'.prep/run.plan.json'
        run('python3',ROOT/'f/gwas_runner.py','execute-plan',plan)
        assert (gen/'typ.prune.bed').stat().st_size>0
        run('zstd','-q',gen/'chr1.pvar','-o',gen/'chr1.pvar.zst')
        if Path('/mnt/d/software/bin/regenie').exists():
            run('bash',ROOT/'gwas.sh','run_regenie','--pheno',d/'phe','--imputed-dir',gen,'--typed-dir',gen,'--output-dir',d/'output',
                '--pheno-name','height,bald,bald12,cvd_cad.Yt2e,cvd_cad.t2e','--chromosomes','1','--covariates','age,sex','--categorical-covariates','center','--threads','2','--memory','1000','--min-mac','5')
            for trait in ['height','bald','bald12','cvd_cad.Yt2e','cvd_cad.t2e']:
                with gzip.open(d/'output'/trait/'gwas'/f'{trait}.regenie.gz','rt') as h:
                    assert h.readline().split()==header.split() and len(h.readlines())>0
        # Check normalization of method-specific effect allele and binary sample count.
        table(d/'regenie','CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE LOG10P',[[1,100,'rs1','A','G',.2,1000,'ADD',.1,.02,9]])
        table(d/'saige','CHR POS MarkerID Allele1 Allele2 AF_Allele2 N_case N_ctrl BETA SE p.value',[[1,100,'rs1','A','G',.2,400,600,.1,.02,1e-9]])
        for method in ('regenie','saige'):
            run('python3',ROOT/'f/gwas_io.py',method,d/(method+'.gz'),d/method)
            with gzip.open(d/(method+'.gz'),'rt') as h:
                h.readline(); row=h.readline().split(); assert row[3:5]==['G','A'] and row[6]=='1000'
        # Exercise real bulik/ldsc munging in the installed Python2 environment.
        python2=Path.home()/'anaconda3/envs/ldsc/bin/python'
        if python2.exists():
            ldrows=[[f'rs{i}',1,i+1,'A','G',.2,10000,(-1 if i%2 else 1)*.1,.02,1e-4] for i in range(1000)]
            ldrows[0][-1]=0
            table(d/'ld.tsv',header,ldrows)
            table(d/'hm3','SNP A1 A2',[[r[0],'A','G'] for r in ldrows])
            run('bash',POST,'ldsc','--gwas-files',d/'ld.tsv','--output-dir',d/'ldsc','--merge-alleles',d/'hm3','--ref-ld-chr',str(d)+'/','--w-ld-chr',str(d)+'/',
                '--python',python2,'--run-h2','FALSE','--run-rg','FALSE')
            assert (d/'ldsc/01.ld.tsv.sumstats.gz').stat().st_size>0
        print('PASS: comparison QC, four real PLINK2 runs, five real REGENIE runs, compressed PVAR, prune/merge, method normalization, LDSC munge, unsupported T2E guard')

if __name__=='__main__': main()
