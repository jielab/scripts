"""A01/A08 regression tests against production functions and public dispatcher.
Run with PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s validation -v.
All synthetic data and outputs live in /tmp.
"""
from pathlib import Path
from unittest import TestCase, main, mock
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import warnings
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def module(name, file):
	if name in sys.modules:
		return sys.modules[name]
	spec = importlib.util.spec_from_file_location(name, ROOT / 'f' / file)
	value = importlib.util.module_from_spec(spec)
	sys.modules[name] = value
	spec.loader.exec_module(value)
	return value

abm = module('c1_abm', 'c1.abm.py')
shared = module('le8_common', '0.common.py')
c5 = module('le8_c5', 'c5.cellulation.py')

class AuditRegression(TestCase):
	def setUp(self):
		self.tmp = tempfile.TemporaryDirectory(prefix='le8-acceptance-', dir='/tmp')
		self.addCleanup(self.tmp.cleanup)
		self.path = Path(self.tmp.name)
		self.clean = mock.patch.dict(os.environ, {k:v for k,v in os.environ.items() if not k.startswith(('LE8_', 'C4_', 'C3_', 'PGS_', 'DATE_FOLLOW_END'))}, clear=True)
		self.clean.start()
		self.addCleanup(self.clean.stop)

	def config(self, *args):
		return abm.s6_validate(abm.configure(abm.reference_parser().parse_args(list(args))))

	def groups(self, n=200):
		p = pd.DataFrame({'eid': [f'{i:05d}' for i in range(n)]})
		m = p.assign(group=[f'{i//2:04d}' for i in range(n)])
		file = self.path / 'families.csv'
		m.to_csv(file, index=False)
		return p, m, file

	def test_environment_family_reaches_outer_and_inner_roles(self):
		p, m, file = self.groups()
		with mock.patch.dict(os.environ, LE8_GROUP_FILE=str(file)):
			a = self.config()
		p = abm.attach_shared_groups(p, a)
		self.assertEqual(p[a.group_col].iloc[0], '0000')
		p['split'] = abm.s6_split(p, a)
		self.assertEqual(p.groupby(a.group_col).split.nunique().max(), 1)
		for part in p.split.unique():
			x = p[p.split.eq(part)].copy()
			x['fold'] = abm.s6_group_folds(abm.s6_groups(x,a), 3, 91)
			self.assertEqual(x.groupby(a.group_col).fold.nunique().max(), 1)

	def test_family_mapping_missing_duplicate_and_conflict_fail(self):
		p, m, file = self.groups()
		a = self.config('--group-file', str(file))
		for bad in [m.iloc[:-1], pd.concat([m, m.iloc[:1]]), m.assign(group=m.group.mask(m.index==0))]:
			bad.to_csv(file,index=False)
			with self.assertRaises(ValueError): abm.attach_shared_groups(p,a)
		m.to_csv(file,index=False)
		a = self.config('--group-file',str(file),'--group-col','family')
		with self.assertRaises(ValueError): abm.attach_shared_groups(p.assign(family='conflicting'),a)

	def test_prepare_joins_group_file_using_string_ids(self):
		p,m,file = self.groups(40)
		p = p.assign(family=m.group, age=55, sex=0, tdi=0, PC1=0, PC2=0, center='A',date_attend='2010-01-01',date_death=None,date_lost='2022-01-01',fod_icd10_cvd_cad=None)
		p.to_csv(self.path/'phenotype.csv',index=False)
		omics=pd.DataFrame({'eid':m.eid,'F1':np.arange(40),'F2':np.arange(40)+1,'F3':np.arange(40)+2})
		omics.to_csv(self.path/'omics.csv',index=False)
		a=self.config('--phe-file',str(self.path/'phenotype.csv'),'--omics-file',str(self.path/'omics.csv'),'--group-file',str(file),'--group-col','family','--residualize','')
		out=self.path/'prepared';out.mkdir()
		actual=abm.prepare(a,out)
		self.assertEqual(actual.eid.iloc[0],'00000')
		self.assertEqual(actual[a.group_col].iloc[0],'0000')

	def test_cutoff_precedence_and_cache_signature(self):
		with mock.patch.dict(os.environ,DATE_FOLLOW_END='2022-12-31',LE8_Y_DATE='custom_date'):
			a=self.config('--demo')
			b=self.config('--demo','--end-date','2024-12-31')
		self.assertEqual(a.diagnosis_col,'custom_date')
		self.assertEqual(a.end_date,'2022-12-31')
		self.assertEqual(b.shared_config_sources['end_date']['overridden_environment'],'2022-12-31')
		self.assertNotEqual(abm.runtime_manifest(a)['signature'],abm.runtime_manifest(b)['signature'])
		phen,_=abm.generate_demo(a);phen,_=abm.outcomes(phen,a);phen['split']=abm.s6_split(phen,a)
		before=abm.cohort_manifest(phen,a);phen.loc[0,'age']+=1
		after=abm.cohort_manifest(phen,a)
		self.assertNotEqual(before['covariates_sha256'],after['covariates_sha256'])
		self.assertEqual(before['test_roster_sha256'],after['test_roster_sha256'])

	def test_shared_outer_roster_and_family_rejection(self):
		p,m,file=self.groups()
		roster=m[['eid']].assign(role=np.where(m.index<40,'test','training'))
		rf=self.path/'outer.csv';roster.to_csv(rf,index=False)
		a=self.config('--group-file',str(file),'--outer-roster',str(rf))
		p=abm.attach_shared_groups(p,a);p['split']=abm.s6_split(p,a)
		self.assertEqual(set(p.loc[p.split.eq('test'),'eid']),set(roster.loc[roster.role.eq('test'),'eid']))
		roster.loc[0,'role']='training';roster.to_csv(rf,index=False)
		with self.assertRaises(ValueError): abm.s6_split(p,a)

	def test_dispatch_only_forwards_c5_options_to_c5(self):
		calls=[]
		with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append((list(map(str,cmd)),env))):
			shared.dispatch_main(['c1_abm,c2_cause,c5_cellulation','--Y','cad','--biom','prot','--analysis-root',str(self.path),
				'--seed','47','--end-date','2024-02-01','--Y-date','custom','--group-file','/tmp/groups.tsv',
				'--contrasts','/tmp/contrast.csv','--matched-draws','100','--cigma-results','/tmp/results.csv',
				'--cigma-cells','/tmp/cells.csv','--allow-untested-cigma','--dry-run'])
		self.assertEqual(len(calls),4)
		for cmd,env in calls:
			self.assertEqual(env['DATE_FOLLOW_END'],'2024-02-01')
			if any('c5.cellulation.py' in x for x in cmd):
				self.assertEqual(cmd[cmd.index('--seed')+1],'47')
				self.assertIn('--contrasts',cmd);self.assertIn('--cigma-cells',cmd)
			else: self.assertNotIn('--contrasts',cmd);self.assertNotIn('--cigma-cells',cmd)
			if any('c1.abm.py' in x for x in cmd): self.assertIn('--group-file',cmd);self.assertIn('--diagnosis-col',cmd)

	def test_abm_dependencies_checked_before_native_scans(self):
		calls=[]
		with mock.patch.dict(os.environ,ABM_PYTHON='/tmp/selected-python'), mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
			shared.dispatch_main(['c1_correlate,c1_abm','--Y','cad,ra','--biom','prot,met','--abm-backend','both',
				'--analysis-root',str(self.path),'--r-bin','/tmp/selected-Rscript','--abm-args','--tree hist','--dry-run'])
		self.assertEqual(len(calls),17)
		for cmd in calls[:8]:
			self.assertEqual(cmd[0],'/tmp/selected-python')
			self.assertIn('--preflight',cmd)
			self.assertIn('/tmp/selected-Rscript',cmd)
			self.assertIn('hist',cmd)
		self.assertIn('c1_correlate',calls[8])
		for check,run in zip(calls[:8],calls[9:]): self.assertEqual(check,run[run.index('--')+1:]+['--preflight'])

	def test_public_abm_defaults_to_cuda_transformer_with_explicit_overrides(self):
		for backend,extra in [('reference',''),('both',''),('reference','--abm-design selective --device cpu')]:
			calls=[]
			with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
				shared.dispatch_main(['c1_abm','--Y','cad','--biom','prot','--analysis-root',str(self.path),
					'--abm-backend',backend,'--abm-args',extra,'--dry-run'])
			for cmd in calls:
				kind=cmd[cmd.index('--backend')+1]
				a=abm.reference_parser().parse_args(cmd[cmd.index('--backend')+2:])
				self.assertEqual(a.device,'cpu' if extra else 'cuda')
				self.assertEqual(a.abm_design,'selective' if extra or kind=='tabicl' else 'selective_attention')

	def test_abm_management_commands_execute_once(self):
		for action in ['evaluate','project','--check-device','--download-model']:
			calls=[]
			with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
				shared.dispatch_main(['c1_abm','--Y','cad','--biom','met','--analysis-root',str(self.path),'--abm-args='+action])
			self.assertEqual(len(calls),1)
			self.assertNotIn('--preflight',calls[0])

	def test_missing_abm_dependency_stops_before_any_native_analysis(self):
		calls=[]
		def fail(cmd,env,dry):
			calls.append(list(map(str,cmd)))
			raise subprocess.CalledProcessError(1,cmd)
		with mock.patch.object(shared,'call',fail),self.assertRaises(subprocess.CalledProcessError):
			shared.dispatch_main(['c1_correlate,c1_abm','--Y','cad','--biom','met','--analysis-root',str(self.path)])
		self.assertEqual(len(calls),1)
		self.assertIn('--preflight',calls[0])
		self.assertFalse(any('0.engine.sh' in c for c in calls[0]))

	def test_explicit_preflight_and_skip_abm(self):
		for flags,n in [(['--preflight'],2),(['--skip-abm','--dry-run'],1)]:
			calls=[]
			with mock.patch.object(shared,'call',lambda cmd,env,dry: calls.append(list(map(str,cmd)))):
				shared.dispatch_main(['c1_correlate,c1_abm','--Y','cad','--biom','met','--analysis-root',str(self.path),*flags])
			self.assertEqual(len(calls),n)
			if '--preflight' in flags: self.assertTrue(all('--preflight' in cmd for cmd in calls))
			else: self.assertIn('c1_correlate',calls[0])

	def test_reader_import_failure_names_selected_interpreter(self):
		phe=self.path/'phenotype.rds';phe.touch()
		omics=self.path/'omics.rds';omics.touch()
		a=self.config('--phe-file',str(phe),'--omics-file',str(omics),'--tree','hist')
		original=abm.importlib.import_module
		for error in [ModuleNotFoundError('pyreadr'),OSError('shared library unavailable')]:
			def broken(name):
				if name=='pyreadr': raise error
				return original(name)
			with mock.patch.object(abm.importlib,'import_module',broken),self.assertRaises(RuntimeError) as caught:
				abm.preflight_inputs(a)
			self.assertIn(sys.executable,str(caught.exception))
			self.assertIn('-m pip install pyreadr',str(caught.exception))

	def test_preflight_rejects_missing_shared_group_and_roster(self):
		for flag in ['--group-file','--outer-roster']:
			phe=self.path/'phenotype.csv';phe.touch()
			omics=self.path/'omics.csv';omics.touch()
			a=self.config('--phe-file',str(phe),'--omics-file',str(omics),flag,str(self.path/'missing.csv'))
			with self.assertRaises(FileNotFoundError): abm.preflight_inputs(a)

	def test_rds_reader_preserves_strings_dates_and_binary_numbers(self):
		file=self.path/'reader.rds'
		script='x <- data.frame(eid=c("00001","00002"), date=as.Date(c("2020-02-29","2021-12-31")), value=c(5e-324,0.125), unused=1:2); saveRDS(x,commandArgs(TRUE)[1])'
		subprocess.run(['Rscript','-e',script,str(file)],check=True)
		for r_bin in ['Rscript','/missing/Rscript']:
			x=abm.read_table(file,columns=['eid','date','value'],r_bin=r_bin)
			self.assertEqual(x.eid.tolist(),['00001','00002'])
			self.assertEqual(pd.to_datetime(x.date).dt.strftime('%Y-%m-%d').tolist(),['2020-02-29','2021-12-31'])
			self.assertEqual(x.value.tolist(),[5e-324,0.125])
			self.assertEqual(x.columns.tolist(),['eid','date','value'])

	def test_abm_signature_survives_transaction_path_changes(self):
		signatures=[]
		published=str(self.path/'published')
		for work in [self.path/'first/results',self.path/'second/results']:
			with mock.patch.dict(os.environ,LE8_ANALYSIS_ROOT=str(work),LE8_PUBLISHED_ROOT=published,LE8_TABLE_WORKSPACE='1'):
				a=self.config('--demo','--analysis-root',str(work),'--run-dir',str(work/'custom'))
				manifest=abm.runtime_manifest(a)
				self.assertEqual(manifest['config']['analysis_root'],published)
				self.assertEqual(manifest['config']['run_dir'],published+'/custom')
				signatures.append(manifest['signature'])
				changed=self.config('--demo','--analysis-root',str(work),'--run-dir',str(work/'custom'),'--end-date','2024-01-01')
				self.assertNotEqual(manifest['signature'],abm.runtime_manifest(changed)['signature'])
		self.assertEqual(signatures[0],signatures[1])

	def test_rds_missing_payload_is_quiet_without_changing_finite_values(self):
		file=self.path/'missing.rds'
		script='x <- data.frame(eid=c("a","b","c"), tdi=c(1,NA_real_,3), category=c(5e-324,1,2)); saveRDS(x,commandArgs(TRUE)[1])'
		subprocess.run(['Rscript','-e',script,str(file)],check=True)
		for r_bin in ['Rscript','/missing/Rscript']:
			x=abm.read_table(file,r_bin=r_bin)
			with warnings.catch_warnings():
				warnings.simplefilter('error',RuntimeWarning)
				self.assertTrue(np.isnan(np.sum(x.tdi.to_numpy())))
				design=abm.S6Clinical(['tdi'],[]).fit(x)
				self.assertTrue(np.isfinite(design.transform(x)).all())
			self.assertEqual(x.category.iloc[0],5e-324)
			np.testing.assert_allclose(x.tdi,[1.,np.nan,3.],equal_nan=True)

	def test_both_preprocessors_preserve_valid_negative_log1p_values(self):
		x=np.array([[-0.073459,0.,1.,np.nan,np.inf]],dtype=np.float32)
		expected=np.log1p(np.where(np.isfinite(x),x,np.nan))
		for cls in [abm.S6Preprocessor,abm.MolecularPreprocessor]:
			prep=cls(.2,[],[],"log1p")
			actual=prep.scale_transform(x)
			np.testing.assert_allclose(actual,expected,equal_nan=True)
			self.assertLess(actual[0,0],0)
			self.assertEqual(x[0,0],np.float32(-0.073459))

	def test_log1p_domain_boundary_and_untransformed_protein_values(self):
		for cls in [abm.S6Preprocessor,abm.MolecularPreprocessor]:
			for value in [-1.,-2.]:
				with self.assertRaisesRegex(ValueError,'x > -1'):
					cls(.2,[],[],"log1p").scale_transform([[value]])
			x=np.array([[-5.,-1.,-.073459,0.,10.]],dtype=np.float32)
			np.testing.assert_array_equal(cls(.2,[],[],"none").scale_transform(x),x)
			near=np.nextafter(np.float32(-1),np.float32(0))
			self.assertTrue(np.isfinite(cls(.2,[],[],"log1p").scale_transform([[near]])).all())

	def test_wide_metabolite_mapping_has_no_fragmentation_or_value_changes(self):
		data={'eid':['00001','00002','00003']}
		data.update({f'p{i}_i0':np.array([i,i+1,i+2],dtype=float) for i in range(315)})
		data['p10_i1']=[999.,999.,999.]
		frame=pd.DataFrame(data,index=[2,4,8])
		mapping=self.path/'mapping.tsv'
		mapping.write_text('data_field\tmet_name\n'+''.join(f'p{i}\tM{i}\n' for i in range(315))+'p1/p0\tRatio\n')
		with warnings.catch_warnings():
			warnings.simplefilter('error',pd.errors.PerformanceWarning)
			mapped,audit=abm.map_metabolites(frame,mapping,'eid')
		self.assertEqual(mapped.index.tolist(),[2,4,8])
		self.assertEqual(mapped.eid.tolist(),data['eid'])
		for i in range(315):np.testing.assert_array_equal(mapped[f'M{i}'],data[f'p{i}_i0'])
		self.assertTrue(np.isnan(mapped.Ratio.iloc[0]))
		self.assertEqual(audit.loc[audit.feature.eq('Ratio'),'nonfinite_to_missing'].item(),1)

	def test_prepare_reports_domain_failure_with_feature_name(self):
		p=pd.DataFrame(dict(eid=['a','b','c'],age=55,sex=0,tdi=0,PC1=0,PC2=0,center='A',
			date_attend='2010-01-01',date_death=None,date_lost=None,fod_icd10_cvd_cad=None))
		p.to_csv(self.path/'phenotype.csv',index=False)
		omics=pd.DataFrame(dict(eid=p.eid,F1=[-.073459,1.,2.],F2=[1.,2.,3.],F3=[2.,3.,4.]))
		omics.to_csv(self.path/'omics.csv',index=False)
		a=self.config('--phe-file',str(self.path/'phenotype.csv'),'--omics-file',str(self.path/'omics.csv'),
			'--residualize','','--transform','log1p')
		out=self.path/'prepared';out.mkdir()
		abm.prepare(a,out)
		self.assertEqual(json.loads((out/'input_audit.json').read_text())['transform_audit']['negative_features'],{'F1':1})
		self.assertLess(np.load(out/'raw.npy')[0,0],0)
		omics.loc[0,'F1']=-1.;omics.to_csv(self.path/'omics.csv',index=False)
		with self.assertRaisesRegex(ValueError,'F1'):abm.prepare(a,out)

	def test_metabolite_negative_qc_preserves_npx_and_valid_observations(self):
		raw=np.array([[-2.,-1.,-.073459,0.,3.,np.nan]],dtype=np.float32)
		original=raw.copy()
		actual=abm.molecular_input_qc(raw,'met')
		np.testing.assert_allclose(actual,[[np.nan,np.nan,np.nan,0.,3.,np.nan]],equal_nan=True)
		np.testing.assert_array_equal(raw,original)
		np.testing.assert_array_equal(abm.molecular_input_qc(raw,'prot'),original)
		np.testing.assert_array_equal(abm.molecular_input_qc(raw,'met',demo=True),original)
		for cls in [abm.S6Preprocessor,abm.MolecularPreprocessor]:
			transformed=cls(.2,[],[],'log1p').scale_transform(actual)
			np.testing.assert_allclose(transformed,np.log1p(actual),equal_nan=True)

	def test_prepare_metabolites_masks_negatives_before_transform(self):
		p=pd.DataFrame(dict(eid=['a','b','c'],age=55,sex=0,tdi=0,PC1=0,PC2=0,center='A',
			date_attend='2010-01-01',date_death=None,date_lost=None,fod_icd10_cvd_cad=None))
		p.to_csv(self.path/'phenotype.csv',index=False)
		omics=pd.DataFrame(dict(eid=p.eid,F1=[-.073459,-1.,-2.],F2=[0.,2.,3.],F3=[2.,3.,4.]))
		omics.to_csv(self.path/'omics.csv',index=False)
		out=self.path/'prepared';out.mkdir()
		for transform in ['none','log1p']:
			a=self.config('--biom','met','--met-input','named','--phe-file',str(self.path/'phenotype.csv'),
				'--omics-file',str(self.path/'omics.csv'),'--residualize','','--transform',transform)
			abm.prepare(a,out)
			audit=json.loads((out/'input_audit.json').read_text())['transform_audit']
			self.assertEqual(audit['negative_features'],{'F1':3})
			self.assertEqual(audit['negative_to_missing_values'],3)
			self.assertEqual(audit['log1p_domain_invalid_values'],0)
			raw=np.load(out/'raw.npy')
			self.assertTrue(np.isnan(raw[:,0]).all())
			np.testing.assert_array_equal(raw[:,1:],omics[['F2','F3']].to_numpy())

	def test_public_shell_seed_environment_and_help(self):
		env=dict(os.environ,PYTHON_BIN=sys.executable,LE8_REPORT_PYTHON=sys.executable,SEED='63')
		r=subprocess.run(['bash',str(ROOT/'le8.sh'),'c5_cellulation','--Y','cad','--biom','prot','--dry-run','--matched-draws','100'],env=env,text=True,capture_output=True,check=True)
		self.assertIn('--seed 63',r.stdout);self.assertIn('--matched-draws 100',r.stdout)

	def test_invalid_c5_settings_fail_at_dispatch(self):
		for args in [['--matched-draws','1'],['--cigma-cells','/tmp/x'],['--universe','/tmp/x']]:
			with self.assertRaises(SystemExit): shared.dispatch_main(['c5_cellulation','--dry-run',*args])

	def cell_inputs(self):
		u=pd.DataFrame({'assay':[f'P{i}' for i in range(12)],'gene':[f'GENE{i}' for i in range(12)],'match_stratum':'same'})
		a=pd.DataFrame({'gene':[f'GENE{i}' for i in range(10)],'cell_type':['A']*5+['B']*5,'source':'fixture','source_version':'v1'})
		p=pd.DataFrame({'model':['YS']*3+['NS']*3,'feature':['P0','P1','P2','P5','P6','P7']})
		c=pd.DataFrame({'model':['YS'],'reference':['NS']})
		for name,d in [('universe',u),('atlas',a),('panels',p),('contrasts',c)]: d.to_csv(self.path/f'{name}.csv',index=False)
		return [item for name in ['universe','atlas','panels','contrasts'] for item in ['--'+name,str(self.path/f'{name}.csv')]]

	def test_c5_dispatch_and_direct_signatures_identical_and_seed_changes(self):
		args=self.cell_inputs();out=self.path/'run'
		base=['--Y','cad','--biom','prot','--analysis-root',str(out),'--no-plots','--seed','57','--matched-draws','100',*args]
		self.assertEqual(c5.main_driver(base),0)
		file=out/'cad/prot/c5_cellulation/c5.completed.json'
		first=json.loads(file.read_text())['signature']
		def execute(cmd,env,dry):
			with mock.patch.dict(os.environ,env,clear=True): self.assertEqual(c5.main_driver(list(map(str,cmd[2:]))),0)
		with mock.patch.object(shared,'call',execute): shared.dispatch_main(['c5_cellulation',*base])
		self.assertEqual(first,json.loads(file.read_text())['signature'])
		base[base.index('--seed')+1]='58';self.assertEqual(c5.main_driver(base),0)
		self.assertNotEqual(first,json.loads(file.read_text())['signature'])

	def test_c5_staged_and_published_evidence_have_identical_provenance(self):
		self.cell_inputs()
		published=self.path/'published';work=self.path/'staging'
		d=pd.DataFrame({'feature':['P0'],'locus':['locus1'],'PP.H4':[.8],'locus_class':['cis']})
		for root in [published,work]:
			f=root/'cad/prot/c3_coloc/c3.coloc_summary.csv';f.parent.mkdir(parents=True);d.to_csv(f,index=False)
		with mock.patch.dict(os.environ,LE8_ANALYSIS_ROOT=str(work),LE8_PUBLISHED_ROOT=str(published)):
			a=c5.stamp(work/'cad/prot/c3_coloc/c3.coloc_summary.csv')
			c5.evidence_long(work/'cad',self.path/'universe.csv',self.path/'staged-evidence')
		b=c5.stamp(published/'cad/prot/c3_coloc/c3.coloc_summary.csv')
		c5.evidence_long(published/'cad',self.path/'universe.csv',self.path/'published-evidence')
		self.assertEqual(a,b)
		for file in ['c5.evidence_long.csv','c5.evidence_sources.csv']:
			self.assertEqual((self.path/'staged-evidence'/file).read_bytes(),(self.path/'published-evidence'/file).read_bytes())

	def test_c5_not_requested_is_explicit(self):
		args=self.cell_inputs();c5.annotate(self.path/'universe.csv',self.path/'atlas.csv',self.path/'panels.csv',self.path/'annotation',plots=False)
		d=pd.read_csv(self.path/'annotation/cell.contrast_status.csv')
		self.assertEqual(d.status.iloc[0],'not_requested')

	def test_c5_cross_fold_contrast_rejected(self):
		self.cell_inputs();p=pd.read_csv(self.path/'panels.csv');p['fold']=['a']*3+['b']*3;p.to_csv(self.path/'panels.csv',index=False)
		pd.DataFrame({'model':['YS|fold=a'],'reference':['NS|fold=b']}).to_csv(self.path/'contrasts.csv',index=False)
		with self.assertRaisesRegex(ValueError,'same fold'):
			c5.annotate(self.path/'universe.csv',self.path/'atlas.csv',self.path/'panels.csv',self.path/'annotation',contrasts=self.path/'contrasts.csv',plots=False)

	def test_r_python_outcome_clock_parity(self):
		a=self.config('--end-date','2020-12-31')
		p=pd.DataFrame(dict(eid=[f'p{i}' for i in range(8)],birth_date='1960-01-01',date_attend=['2020-01-01']*7+['2020-12-31'],
			date_death=[None,None,None,None,'2020-06-01',None,None,None],date_lost=[None]*5+['2020-04-01',None,None],
			fod_icd10_cvd_cad=['2019-12-31','2020-01-01','2020-12-31','2021-01-01','2020-07-01','2020-04-01',None,None]))
		file=self.path/'outcomes.csv';p.to_csv(file,index=False)
		r=subprocess.run(['Rscript',str(ROOT/'validation/audit_acceptance.R'),'--outcomes',str(file)],env=dict(os.environ,DATE_FOLLOW_END=a.end_date),text=True,capture_output=True,check=True)
		data=json.loads(r.stdout);out,audit=abm.outcomes(p,a)
		self.assertEqual(data['prevalent'],out.prevalent.astype(int).tolist())
		self.assertEqual(data['event'],[int(x) if valid else None for x,valid in zip(out.event,out.endpoint_valid)])
		np.testing.assert_allclose([np.nan if x is None else x for x in data['time']],out.time.where(out.endpoint_valid).to_numpy(),equal_nan=True)

if __name__=='__main__': main()
