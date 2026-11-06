"""REGENIE X metadata and generated-command regressions, without user data."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import test_runtime

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('audit', ROOT/'f/gwas.audit.py')
audit=importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)

class XAudit(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup); self.p=Path(self.tmp.name)
  for name,content in {'phe':'FID IID bald\na a 1\nb b 2\n', 'fam':'a a 0 0 1 -9\nb b 0 0 2 -9\n', 'psam':'#FID IID SEX\na a 1\nb b 2\n', 'pvar':'#CHROM POS ID REF ALT\nPAR1 100000 a A G\nX 5000000 b A G\nPAR2 155000000 c A G\n'}.items():
   (self.p/name).write_text(content)
 def run_audit(self,**kwargs):
  return audit.audit(*(self.p/n for n in ['phe','fam','psam','pvar']),**kwargs)
 def test_male_subset_and_regions(self):
  result=self.run_audit(sex='male'); self.assertEqual(result['audited_overlap'],1)
  self.assertEqual(result['regions'],dict(PAR1=1,nonPAR=1,PAR2=1))
  self.assertFalse(result['genotype_encoding_checked'])
 def test_unknown_and_discordant_sex(self):
  for sex in ['0','2']:
   (self.p/'psam').write_text('#FID IID SEX\na a '+sex+'\nb b 2\n')
   with self.assertRaises(ValueError): self.run_audit()
 def test_non_x_and_bad_par_rejected(self):
  for row in ['Y 100000','PAR1 5000000','X 160000000']:
   (self.p/'pvar').write_text('#CHROM POS ID REF ALT\n'+row+' rs A G\n')
   with self.assertRaises(ValueError): self.run_audit()
 def test_keep_filter(self):
  (self.p/'keep').write_text('FID IID\nb b\n')
  self.assertEqual(self.run_audit(keep=self.p/'keep')['audited_overlap'],1)
  with self.assertRaises(ValueError): self.run_audit(keep=self.p/'keep',sex='male')

class RegenieCommands(unittest.TestCase):
 def setUp(self):
  self.h=test_runtime.PortableRuntime(); self.h.setUp(); self.addCleanup(self.h.doCleanups)
 def test_trait_commands_and_default_chr(self):
  h=self.h; geno,_=h.genotype_fixture()
  for chrom in ['X','Y']:
   for ext in ['pgen','pvar','psam']: (geno/f'chr{chrom}.{ext}').write_text('fixture\n')
  for kind,trait in [('ordinal','bald'),('bt','bald12'),('t2e','cvd_cad.t2e')]:
   out=h.work/kind
   h.call('gwas.sh','regenie','--submit','local','--run-cmd','FALSE','--dir-imputed',geno,'--dir-typ',geno,'--dir-output',out,'--pheno',h.work/'phe','--pheno-name',trait,'--type',kind,'--sex','male')
   plan=json.loads((out/trait/'gwas/regenie/run.plan.json').read_text())
   jobs={j['name']:j for j in plan['jobs']}
   self.assertEqual(len([k for k in jobs if k.endswith('.association')]),23)
   self.assertIn('genotype-audit',jobs); self.assertIn('chrX.association',jobs); self.assertNotIn('chrY.association',jobs)
   self.assertEqual(plan['config']['covariates'],'age,sex,PC1,PC2,PC3,PC4')
   for name in ['null-model','chrX.association']:
    cmd=jobs[name]['cmd']; self.assertEqual(cmd[cmd.index('--sex-specific')+1],'male')
    self.assertEqual(cmd[cmd.index('--par-region')+1],'b37'); self.assertIn('--strict',cmd)
    if kind=='ordinal': self.assertIn('--force-qt',cmd)
    elif kind=='bt': self.assertIn('--bt',cmd)
    else:
     self.assertIn('--t2e',cmd); self.assertEqual(cmd[cmd.index('--eventColList')+1],'cvd_cad.Yt2e')
 def test_incomplete_sync_fails_before_tools(self):
  h=self.h; geno,plink=h.genotype_fixture()
  (geno/'chr1.pgen').rename(geno/'chr1.pgen.filepart')
  (geno/'chr22.pvar').unlink()
  result=h.call('gwas.sh','regenie','--submit','local','--dir-imputed',geno,'--dir-typ',geno,'--dir-output',h.work/'out','--pheno',h.work/'missing.phe','--pheno-name','bald','--type','ordinal','--plink2',plink,success=False)
  self.assertIn('INPUT_TRANSFER_IN_PROGRESS:',result.stderr)
  self.assertIn(str(geno/'chr22.pvar'),result.stderr)
  self.assertIn(str(geno/'chrX.pgen'),result.stderr)
  self.assertIn('No GWAS tools were started',result.stderr)
  self.assertFalse(list(geno.glob('*.calls')))
 def test_separate_covariates(self):
  h=self.h; geno,_=h.genotype_fixture(); out=h.work/'cimt'
  h.call('gwas.sh','regenie','--submit','local','--run-cmd','FALSE','--dir-imputed',geno,'--dir-output',out,'--pheno',h.work/'cimt.phe','--pheno-name','phe1,phe2','--type','qt','--continuous-covariate-file',h.work/'cimt.cont.cov','--categorical-covariate-file',h.work/'cimt.cat.cov')
  for trait in ['phe1','phe2']:
   cfg=json.loads((out/trait/'gwas/regenie/run.plan.json').read_text())['config']
   self.assertEqual(cfg['covariates'],'*'); self.assertEqual(cfg['categorical_covariates'],'*')
if __name__=='__main__': unittest.main()
