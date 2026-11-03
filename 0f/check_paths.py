#!/usr/bin/env python3
"""Check selected projects' default input inventory, without reading large files.

This inventory covers default source roots and key references. Each scientific
module additionally validates its selected traits, samples and complete inputs.
"""
import argparse
import json
import os
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project',choices=['all','allele3','grid','gu','gwas','ld','le8','metag'],default='all')
p.add_argument('--chr',default='1-22',help='Autosomes for genotype inventory, comma/range syntax')
p.add_argument('--json',action='store_true')
a=p.parse_args()
chroms=[]
for token in a.chr.split(','):
 lo,hi=map(int,token.split('-')) if '-' in token else (int(token),int(token))
 if not 1<=lo<=hi<=22:p.error('--chr must be within 1..22')
 chroms+=list(range(lo,hi+1))
roots={k:Path(os.environ[k]) for k in ('DIR_DAT','DIR_BIG','DIR_TMP','DIR_GEN','DIR_ANNOT','DIR_GWAS','DIR_LD','DIR_UKB','DIR_FILES','DIR_METAG')}
rows=[]
def check(project,*alternatives):
 if a.project not in ('all',project):return
 paths=list(map(Path,alternatives));ok=any(x.is_dir() or (x.is_file() and x.stat().st_size>0) for x in paths)
 rows.append(dict(project=project,present=ok,expected=[str(x) for x in paths]))
for c in dict.fromkeys(chroms):
 for project,dataset,build in [('ld','1kg','38'),('gu','1kg','37')]:
  prefix=roots['DIR_GEN']/dataset/build/'pfile'/f'chr{c}'
  for ext in ('.pgen','.psam'):check(project,str(prefix)+ext)
  check(project,str(prefix)+'.pvar.zst',str(prefix)+'.pvar')
check('ld',roots['DIR_GEN']/'1kg/38/samples.txt')
check('gu',roots['DIR_GEN']/'archaic/37/vcf')
check('gu',roots['DIR_ANNOT']/'axt/37')
check('gu',roots['DIR_GEN']/'1hgp/37/GRCH37.fasta.gz')
check('allele3',roots['DIR_ANNOT']/'gnomad/sites')
check('allele3',roots['DIR_GEN']/'1hgp/38/GRCH38.fasta.gz')
check('grid',roots['DIR_UKB']/'Rdata/phe.rds')
check('grid',roots['DIR_FILES']/'DiscoDivas/med.g1000.4pop.tsv')
check('grid',roots['DIR_FILES']/'DiscoDivas/g1k_hm3_maf5_woamb_wolr.pca.weight')
check('grid',roots['DIR_LD']/'csx/snpinfo_mult_1kg_hm3')
check('grid',roots['DIR_GWAS']/'4grid/common')
check('gwas',roots['DIR_UKB']/'common/ukb.phe')
check('gwas',roots['DIR_UKB']/'common/snp.lst')
check('le8',roots['DIR_UKB']/'Rdata/phe.rds')
check('le8',roots['DIR_GWAS']/'main')
check('le8',roots['DIR_GWAS']/'prot')
check('metag',Path(os.environ['METAG_REFERENCE_RAW']))
check('metag',Path(os.environ['METAG_READS_RAW']))
if a.json:print(json.dumps({'roots':{k:str(v) for k,v in roots.items()},'checks':rows},indent=2))
else:
 for r in rows: print(('present ' if r['present'] else 'missing ')+ ' OR '.join(r['expected']) + ' ['+r['project']+']')
 print('Inventory covers default roots; run the selected module preflight for trait/sample-specific checks.')
raise SystemExit(2 if any(not x['present'] for x in rows) else 0)
