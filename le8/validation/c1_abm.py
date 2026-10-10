"""C1 ABM: preprocessing, retrieval, training, devices, simulations and public workflow.

Run: python -m pytest validation/c1_abm.py
Fixtures: python validation/c1_abm.py --smoke-data --out /tmp/le8-inputs
"""
import argparse
import copy
import csv
import dataclasses
import hashlib
import importlib.util
import json
import numpy as np
import os
import pandas as pd
import pickle
import pytest
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import types
import warnings
from dataclasses import replace
from pathlib import Path
from threadpoolctl import threadpool_limits
from unittest import TestCase, mock
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

def load_module(name, filename):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / 'f' / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

m = abm = load_module('le8_validation_abm', 'c1.abm.py')
shared = load_module('le8_validation_dispatch', '0.common.py')
r = m.s7_resources()
BlockedKNN = r.BlockedKNN


def config(*extra):
	return m.s7_validate(m.configure(m.reference_parser().parse_args([
		'--demo','--abm-design','selective_attention','--tree','hist','--tree-estimators','3',
		'--epochs','1','--pretrain-epochs','0','--width','8','--heads','2','--layers','1','--tokens','3',
		'--batch-size','128','--selective-folds','3','--selective-min-events','2',
		'--s7-radius-ref-n','20','--s7-neighbors-per-view','4','--s7-bandwidth-grid','1',
		'--covariates','age,sex','--categorical','sex','--residualize','','--horizon','5',
		'--device','cpu','--cores','2','--group-col','family','--group-namespace','synthetic',
		'--bootstrap','0','--selective-gate-trees','10','--c-grid','.01',*extra])))

def synthetic(n=2000):
	rng=np.random.default_rng(147)
	x=rng.normal(size=(n,8));x[:,5]=(rng.random(n)<.01)*20
	p=pd.DataFrame(dict(eid=[f'synth-{i:05d}' for i in range(n)],family=[f'family-{i//2:05d}' for i in range(n)],
		age=rng.uniform(40,70,n),sex=rng.integers(0,2,n)))
	p['event']=(rng.random(n)<m.expit(-.5+.5*x[:,0])).astype(int)
	p['time']=np.where(p.event,rng.uniform(.2,4.8,n),rng.uniform(6,9,n))
	return x,p,[f'F{i}' for i in range(8)]

@pytest.fixture(scope='module')
def fitted():
	with tempfile.TemporaryDirectory(prefix='s7-production-',dir='/tmp') as tmp,threadpool_limits(2):
		m.torch.set_num_threads(2);a=config();x,p,names=synthetic();p['role']=m.s7_split(p,a)
		out=Path(tmp)/'fit';b=m.s7_fit(x,p,names,a,out)
		yield dict(a=a,x=x,p=p,names=names,b=b,out=out,root=Path(tmp))

def request_fixture(fitted):
	p=fitted['p'].iloc[:8][['eid','age','sex']].copy();p['eid']=['external-'+str(i) for i in range(len(p))]
	groups=pd.DataFrame(dict(eid=p.eid,group=['external-family-'+str(i//2) for i in range(len(p))]))
	root=fitted['root'];phe=root/'external.csv';gf=root/'groups.csv';omics=root/'omics.csv'
	p.to_csv(phe,index=False);groups.to_csv(gf,index=False)
	z=pd.DataFrame(fitted['x'][:8],columns=fitted['names']);z.insert(0,'eid',p.eid);z.to_csv(omics,index=False)
	request=m.ProjectionRequest(str(phe),str(omics),str(root/'project.csv'),str(gf),group_namespace='synthetic')
	return p,groups,z,request

def test_T01_T07_T48_projection_cli_and_frozen_cpu(fitted):
	p,g,z,r=request_fixture(fitted)
	command=[sys.executable,str(ROOT/'f/c1.abm.py'),'project','--run-dir',str(fitted['out']),
		'--phe-file',r.phe_file,'--omics-file',r.omics_file,'--output',r.output,
		'--projection-group-file',r.group_file,'--projection-group-namespace','synthetic',
		'--device','cpu','--cores','2','--memory-limit-gb','0']
	subprocess.run(command,check=True,capture_output=True,text=True)
	one=pd.read_csv(r.output).set_index('eid')
	assert (one.validation_status=='external_validation').all()
	p.iloc[::-1].to_csv(r.phe_file,index=False);g.sample(frac=1,random_state=5).to_csv(r.group_file,index=False);z.iloc[::-1].to_csv(r.omics_file,index=False)
	two=m.s7_project(fitted['out'],r).set_index('eid').loc[one.index]
	np.testing.assert_allclose(one.policy_risk,two.policy_risk,atol=1e-6,rtol=0)
	learner=fitted['b']['models']['attention_retrieval_weighted'];state=learner.export_state()
	restored=m.S7AttentionLearner.restore_state(state,'cpu');meta,_=m.prepare_projection_metadata(fitted['b']['projection_schema'],r,z.eid.tolist())
	np.testing.assert_allclose(learner.predict(fitted['x'][:8],meta).probability,restored.predict(fitted['x'][:8],meta).probability,atol=1e-7)

@pytest.mark.parametrize('case',['missing','conflicting','participant','family','namespace_missing'])
def test_T02_T03_T04_T05_projection_roster(fitted,case):
	p,g,z,r=request_fixture(fitted)
	if case=='missing':g=g.iloc[:-1]
	if case=='conflicting':p['source']=['wrong']*len(p);r=replace(r,group_col='source')
	if case=='participant':
		dev=fitted['p'].loc[fitted['p'].role.eq('tune_gate')].iloc[0]
		old=p.eid.iloc[0];p.loc[0,'eid']=dev.eid;g.loc[0,'eid']=dev.eid;z.loc[0,'eid']=dev.eid
		r=replace(r,group_namespace='different_namespace') # Namespace cannot bypass participant exclusion.
	if case=='family':g.loc[0,'group']=next(iter(fitted['b']['development_groups']))
	if case=='namespace_missing':r=replace(r,group_namespace=None)
	p.to_csv(r.phe_file,index=False);g.to_csv(r.group_file,index=False);z.to_csv(r.omics_file,index=False)
	with pytest.raises(ValueError):m.load_projection(fitted['b'],r)

def test_T06_string_ids_csv_rds(tmp_path):
	p=pd.DataFrame(dict(eid=['00001','00002'],age=[50.,60.]));csv=tmp_path/'ids.csv';p.to_csv(csv,index=False)
	rds=tmp_path/'ids.rds'
	subprocess.run(['Rscript','-e',"args<-commandArgs(TRUE);saveRDS(data.frame(eid=c('00001','00002'),age=c(50,60)),args[1])",str(rds)],check=True)
	assert m.read_table(csv).eid.tolist()==m.read_table(rds).eid.tolist()==['00001','00002']

def test_T11_T12_technical_tail_and_schema():
	a=config('--residualize','sex');x,p,n=synthetic(1000)
	x[:,5]=0;x[:3,5]=20
	prep=m.S7Preprocessor(vars(a)).fit(x,m.s7_metadata(p,a),m.s7_schema(n,a))
	before=m.s7_state_hash(prep);xx=x[:2].copy();xx[0,0]=1e300;xx[1,0]=np.nan
	z=prep.transform(xx,p.iloc[:2]);j=list(prep.tail.keep_).index(0)
	assert z.tail_cap_hit[0,j] and not z.observed[1,j] and z.tail_up[1,j]==0
	assert before==m.s7_state_hash(prep)
	assert z.design().shape[1]==4*len(z.feature_names)+z.clinical.shape[1]
	p2=p.iloc[:2].copy();p2['sex']=99
	assert prep.transform(xx,p2).unknown_category.all()
	one=m.S7Preprocessor(vars(a)).fit(x[:,[0]],m.s7_metadata(p,a),m.s7_schema(['one_assay'],a))
	assert one.transform(x[:2,[0]],p.iloc[:2]).bulk.shape==(2,1)

def test_T17_T18_geometry_stable_candidates(fitted):
	learner=fitted['b']['models']['attention_retrieval_full'];ix=np.flatnonzero(fitted['p'].role.eq('test'))[:5]
	batch=learner.prep.transform(fitted['x'][ix],fitted['p'].iloc[ix]);z,_=m.s7_encode(learner.network,batch)
	before=m.s7_state_hash(learner.memory.geometry)
	a=learner.memory.geometry.candidates(batch,z)
	b=learner.memory.geometry.candidates(batch.take(np.arange(5)[::-1]),z[::-1])
	np.testing.assert_array_equal(a[0],b[0][::-1]);np.testing.assert_allclose(a[1],b[1][::-1])
	assert before==m.s7_state_hash(learner.memory.geometry)
	assert a[0].shape[1]<=3*fitted['a'].s7_neighbors_per_view
	batch.tail_up.fill(0);batch.tail_down.fill(0)
	c=learner.memory.geometry.candidates(batch,z)
	assert c[2].any()
	# All donors equidistant: a bank reorder must preserve the selected identities.
	cfg={**vars(fitted['a']),'s7_views':['learned']}
	bank=learner.prep.transform(fitted['x'][:100],fitted['p'].iloc[:100])
	emb=np.tile([-1.,1.],(100,1));emb[50:]*=-1
	query=batch.take(np.array([0]));qz=np.zeros((1,2))
	left=m.S7RetrievalGeometry().fit(bank,emb,cfg)
	order=np.random.default_rng(10).permutation(100)
	right=m.S7RetrievalGeometry().fit(bank.take(order),emb[order],cfg)
	left.radii['learned']=right.radii['learned']=1.
	li,_,la,_=left.candidates(query,qz);ri,_,ra,_=right.candidates(query,qz)
	assert la.any() and ra.any()
	np.testing.assert_array_equal(left.ids[li[la]],right.ids[ri[ra]])

def test_T19_T20_masking_and_sdpa(fitted):
	learner=fitted['b']['models']['attention_direct_full'];b=learner.prep.transform(fitted['x'][:3],fitted['p'].iloc[:3])
	data=b.tensors('cpu');data['observed'][:]=False
	learner.network.eval()
	with m.torch.no_grad():
		z1=learner.network(**data)[0]
		data['bulk'][:]=100;data['tail_up'][:]=1000;data['tail_down'][:]=12
		z2=learner.network(**data)[0];z3=learner.network(**data)[0]
	np.testing.assert_array_equal(z1.numpy(),z2.numpy());np.testing.assert_array_equal(z2.numpy(),z3.numpy())

def test_T21_T22_T23_gradients_and_values(fitted):
	model=fitted['b']['models']['attention_retrieval_full']
	assert {'query.weight','key.weight'}<=model.gradient_names
	assert model.fit_status['supervised_steps']>=2
	assert model.fit_status['selected_tensor_changed_count']>0 or model.selected_epoch==0
	ix=np.flatnonzero(fitted['p'].role.eq('test'))[:10];b=model.prep.transform(fitted['x'][ix],fitted['p'].iloc[ix])
	with m.torch.no_grad():
		q,logit,_=model.network(**b.tensors('cpu'));p=m.torch.sigmoid(logit)
		jj,r,allow,_=model.memory.geometry.candidates(b,q.numpy())
		left,dl=m.s7_torch_borrow(model.network,q,p,model.memory,jj,r,allow,1,model.config)
		memory=copy.deepcopy(model.memory);memory.labels=1-memory.labels
		right,dr=m.s7_torch_borrow(model.network,q,p,memory,jj,r,allow,1,model.config)
		np.testing.assert_array_equal(dl['logits'].numpy(),dr['logits'].numpy())
		np.testing.assert_array_equal(dl['weights'].numpy(),dr['weights'].numpy())
		assert np.max(np.abs(left.numpy()-right.numpy()))>1e-7
	assert model.prep.feature_names==fitted['b']['models']['attention_direct_full'].prep.feature_names

def test_T24_T31_batch_and_future_columns(fitted):
	b=fitted['b'];ix=np.flatnonzero(fitted['p'].role.eq('test'))[:100];p=fitted['p'].iloc[ix];x=fitted['x'][ix]
	before=m.s7_freeze_components(b)
	allp=m.s7_predict(b,x,p);one=m.s7_predict(b,x[:1],p.iloc[:1]);rev=m.s7_predict(b,x[::-1],p.iloc[::-1])
	np.testing.assert_allclose(one.policy_risk,allp.policy_risk.iloc[:1],atol=1e-6,rtol=0)
	np.testing.assert_allclose(allp.policy_risk,rev.policy_risk.iloc[::-1],atol=1e-6,rtol=0)
	p=p.copy();p['time']=999.;p['event']=1-p.event;p['future_death']='2099-01-01'
	poison=m.s7_predict(b,x,p)
	np.testing.assert_array_equal(allp.policy_risk,poison.policy_risk)
	assert before==m.s7_freeze_components(b)

@pytest.mark.parametrize('role',['test','tune_gate','calibration_audit'])
def test_T25_T27_T28_role_outcomes_do_not_refit_candidates(fitted,role):
	a=fitted['a'];p=fitted['p'].copy();ix=np.flatnonzero(p.role.eq(role));rng=np.random.default_rng(31)
	p.loc[ix,['time','event']]=p.loc[rng.permutation(ix),['time','event']].to_numpy()
	with threadpool_limits(2):other=m.s7_fit(fitted['x'],p,fitted['names'],a,fitted['root']/('poison-'+role))
	assert fitted['b']['component_hashes']['models']==other['component_hashes']['models']
	assert fitted['b']['component_hashes']['calibrators']==other['component_hashes']['calibrators']
	if role!='tune_gate':
		assert fitted['b']['component_hashes']['release_gate']==other['component_hashes']['release_gate']
		assert fitted['b']['component_hashes']['coverage_rules']==other['component_hashes']['coverage_rules']
	if role=='test':
		test=np.flatnonzero(p.role.eq('test'))
		old=m.s7_predict(fitted['b'],fitted['x'][test],p.iloc[test]);new=m.s7_predict(other,fitted['x'][test],p.iloc[test])
		np.testing.assert_array_equal(old.policy_risk,new.policy_risk)
		np.testing.assert_array_equal(old.released,new.released)

def test_T26_recipient_family_poisoning_full_neural_pipeline(fitted):
	a=fitted['a'];ix=np.flatnonzero(fitted['p'].role.eq('build'));p=fitted['p'].iloc[ix].reset_index(drop=True)
	data=m.s7_fit_data(fitted['x'][ix],p,m.s7_schema(fitted['names'],a),a,'build')
	fold=m.s6_group_folds(m.s6_groups(p,a),a.selective_folds,a.seed+1200);u=np.flatnonzero(fold==0)
	poison=copy.deepcopy(data);poison.event[u]=1-poison.event[u];poison.time[u]=np.where(poison.event[u],2.,8.)
	with threadpool_limits(2):
		left,lv=m.s7_crossfit_training_gate(data,vars(a),only_recipient_fold=0)
		right,rv=m.s7_crossfit_training_gate(poison,vars(a),only_recipient_fold=0)
	np.testing.assert_array_equal(lv['weight'][u],rv['weight'][u]);np.testing.assert_array_equal(lv['hard_mask'][u],rv['hard_mask'][u])
	assert left.audit.pilot_state_hash.tolist()==right.audit.pilot_state_hash.tolist()

def test_T29_T38_T41_T42_evaluation_no_fit(fitted):
	path=fitted['out']/'model_bundle.joblib';before=m.hashlib.sha256(path.read_bytes()).hexdigest()
	one=m.s7_evaluate(fitted['out']);two=m.s7_evaluate(fitted['out'])
	assert before==m.hashlib.sha256(path.read_bytes()).hexdigest()
	pd.testing.assert_frame_equal(one,two)
	assert one.groupby('subset').N.nunique().max()==1
	assert one.groupby('subset').observed_cases.nunique().max()==1
	curve=pd.read_csv(fitted['out']/'coverage_curve.csv')
	assert {'Uno_C_horizon','weighted_pair_denominator','comparable_pairs'}<=set(curve)
	assert curve.loc[curve.requested_coverage.eq(1),'actual_coverage'].eq(1).all()
	pairs=pd.read_csv(fitted['out']/'paired_contrasts.csv')
	assert {'attention_effect','retrieval_effect','training_weight_effect','weight_specificity','combined_method','policy_effect'}==set(pairs.contrast)

def test_T30_public_family_argument(tmp_path):
	cmd=['bash',str(ROOT/'le8.sh'),'c1_abm','--Y','cad','--biom','prot','--group-file','/tmp/groups.csv',
		'--dry-run','--abm-args','--abm-design selective_attention --device cpu']
	r=subprocess.run(cmd,capture_output=True,text=True,check=True)
	assert '--group-file /tmp/groups.csv' in r.stdout and '--abm-design selective_attention' in r.stdout
	a=config('--group-col','','--group-file','/tmp/groups.csv');assert a.group_col=='.le8_family'

def test_T32_module_provenance(tmp_path):
	file=tmp_path/'modules.tsv';pd.DataFrame(dict(feature=['F0'],module=['M0'])).to_csv(file,sep='\t',index=False)
	a=config('--module-file',str(file))
	with pytest.raises(ValueError):m.s7_module_file(a,['synthetic'])
	pd.DataFrame(dict(feature=['F0'],module=['M0'],provenance=['external_fixed'])).to_csv(file,sep='\t',index=False)
	assert m.s7_module_file(a,['synthetic'])==str(file)

def test_T33_T34_T35_T43_release_contract(fitted):
	b=fitted['b'];a=fitted['a'];ix=np.flatnonzero(fitted['p'].role.eq('test'));p=fitted['p'].iloc[ix]
	assert m.S7CoverageRule().fit(np.arange(3),['a','b','c'],1,1).apply(np.array([-1e9,0,1e9]),['x','y','z']).all()
	assert b['release_gate'].target=={'primary':a.s7_primary,'reference':'baseline_ensemble_full'}
	assert b['audit_state']['status']=='insufficient_audit_information'
	assert b['audit_state']['minimum_valid_boot']==200
	# Actual gate fitted to a nonpositive conditional target, no invented release.
	z=pd.DataFrame(np.zeros((100,len(m.S7_GATE_COLUMNS))),columns=m.S7_GATE_COLUMNS)
	y=np.arange(100)%2
	gate=m.S7ReleaseGate().fit(z,y,np.ones(100),.8-.6*y,.1+.8*y,vars(a))
	assert (gate.predict(z)['gain']<0).all()
	out=m.s7_predict(b,fitted['x'][ix],p)
	assert not out.released.any();np.testing.assert_array_equal(out.policy_risk,out[b['fallback']])

def test_T36_random_controls_match_distribution():
	values=dict(clinical_risk=np.linspace(0,1,100),weight=np.linspace(.2,1,100),hard_mask=np.arange(100)%3==0)
	w,h=m.s7_random_controls(values,12)
	np.testing.assert_array_equal(np.sort(w),np.sort(values['weight']));assert h.sum()==values['hard_mask'].sum()
	assert not np.array_equal(w,values['weight'])

def test_T47_numerical_manifest_and_report_only_changes():
	source=(ROOT/'f/c1.abm.py').read_text()
	report=source.replace("title='Frozen S7 models: common technically eligible participants'", "title='Updated figure caption'")
	assert source!=report
	assert m.s7_numerical_code_hash(source)==m.s7_numerical_code_hash(report)
	assert m.s7_numerical_code_hash(source)!=m.s7_numerical_code_hash(source.replace('a.s7_prior_strength < 0','a.s7_prior_strength < 1'))
	a=config();original=m.s7_manifest(a)['signature']
	a.s7_max_borrow=.7
	assert m.s7_manifest(a)['signature']!=original

def test_prespecified_conditional_censoring_is_build_only(fitted):
	b=fitted['b'];ix=np.flatnonzero(fitted['p'].role.eq('build'));a=fitted['a']
	p=fitted['p'].iloc[ix]
	independent=m.S7ConditionalCensoring().fit(m.s7_metadata(p,a),p[['time','event']],vars(a))
	assert m.s7_state_hash(independent)==m.s7_state_hash(b['conditional_censoring'])
	assert independent.columns==['age','sex']

def test_T37_filtered_donor_bank_excludes_labels(fitted):
	a=fitted['a'];x,p,n=synthetic(600);schema=m.s7_schema(n,a);data=m.s7_fit_data(x,p,schema,a,'build')
	mask=np.ones(450);mask[::3]=0
	with threadpool_limits(2):model=m.S7AttentionLearner().fit(data.take(np.arange(450)),data.take(np.arange(450,600)),mask,{**vars(a),'supervised_filter':True})
	assert not set(p.eid.iloc[np.flatnonzero(mask==0)]) & set(model.memory.ids)

def test_T39_uno_naive_pairs_and_ties():
	p=pd.DataFrame(dict(time=[1.,2.,2.,3.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.],event=[1,1,0,1]+[0]*10))
	km=m.S6Censoring().fit(p,2.5,.05);risk=np.array([.4,.4,.7,.1]+[.2]*10)
	num=den=0.;count=0
	for i in range(len(p)):
		for j in range(len(p)):
			if p.event.iloc[i] and p.time.iloc[i]<=km.horizon and p.time.iloc[j]>p.time.iloc[i]:
				w=1/km.at([p.time.iloc[i]],left=True)[0]**2;den+=w;count+=1;num+=w*(int(risk[i]>risk[j])+.5*int(risk[i]==risk[j]))
	actual=m.s7_uno(p.time.to_numpy(),p.event.to_numpy(),risk,km)
	assert abs(actual['Uno_C_horizon']-num/den)<1e-10
	assert actual['comparable_pairs']==count and abs(actual['weighted_pair_denominator']-den)<1e-10

def test_T40_paired_bootstrap_uses_same_family_draw(fitted):
	p=fitted['p'].iloc[:100];prob=np.linspace(.1,.8,100)
	tab=m.s7_paired_bootstrap(p.time.to_numpy(),p.event.to_numpy(),{'a':prob,'b':prob},np.ones(100,bool),p.family.to_numpy(),fitted['b']['censoring'],[('a','b')],30,1)
	assert tab.delta.eq(0).all() and tab.lower.eq(0).all() and tab.upper.eq(0).all()
	assert set(tab.metric)=={'AUC_IPCW','Uno_C_horizon','Brier_IPCW','LogLoss_IPCW'}

def test_T44_T46_mode_paths_and_signature():
	a=config();b=config('--seed','7');c=config('--horizon','6')
	assert len({m.runtime_manifest(v)['signature'] for v in (a,b,c)})==3
	assert m.output_directory(a).name=='abm_selective_attention'
	default=m.configure(m.reference_parser().parse_args([]))
	assert default.abm_design=='selective_attention' and default.device=='cuda'
	old=m.configure(m.reference_parser().parse_args(['--abm-design','selective','--device','cpu']))
	assert old.abm_design=='selective' and m.output_directory(old).name=='abm_reference'

@pytest.mark.parametrize('backend',['tabicl','both'])
def test_T45_T52_backend_rejected_before_training(backend):
	cmd=['bash',str(ROOT/'le8.sh'),'c1_abm','--abm-backend',backend,'--dry-run','--abm-args','--abm-design selective_attention']
	r=subprocess.run(cmd,capture_output=True,text=True)
	assert r.returncode!=0 and 'incompatible' in r.stderr+r.stdout

def test_T51_explicit_unavailable_cuda_fails():
	# Hide the physical GPU in a real preflight subprocess; never substitute a fake network.
	cmd=[sys.executable,str(ROOT/'f/c1.abm.py'),'--abm-design','selective_attention','--device','cuda','--check-device','--memory-limit-gb','0']
	r=subprocess.run(cmd,env={**os.environ,'CUDA_VISIBLE_DEVICES':''},capture_output=True,text=True)
	assert r.returncode!=0 and 'CUDA requested but unavailable' in r.stderr+r.stdout

def test_resume_preserves_exact_numeric_contract_across_dispatcher_changes():
	source=(ROOT/'f/c1.abm.py').read_text()
	contract=m.s7_resume_code_hash(source)
	assert contract=='c8f60c3a5b9d8c3c973ba9918fd98897906ad8c140f6556239ba7a7180408134'
	old_hint = 'LightGBM missing; install requirements.txt or explicitly choose --tree hist before training'
	new_hint = 'LightGBM missing; run ./install.sh --abm or explicitly choose --tree hist before training'
	assert m.s7_resume_code_hash(source.replace(new_hint, old_hint)) == contract
	old=dict(signature='old',config=dict(seed=2026),inputs=[dict(sha256='input')],
		code=dict(abm_numeric_sha256='4e23ce5e9ef7baefca374fb51acf0a920cf0692551341ed60e038e4d0d1d3a99',
			retrieval_kernel_sha256='kernel',dispatcher=['old dispatcher']))
	new=copy.deepcopy(old);new['signature']='new'
	new['code'].update(abm_numeric_sha256='new source',abm_resume_sha256=contract,dispatcher=['new dispatcher'])
	assert m.s7_manifests_compatible(old,new)
	for path,value in [(('config','seed'),2027),(('code','retrieval_kernel_sha256'),'changed kernel'),
		(('code','abm_resume_sha256'),'changed model')]:
		changed=copy.deepcopy(new);changed[path[0]][path[1]]=value
		assert not m.s7_manifests_compatible(old,changed)
	changed=copy.deepcopy(new);changed['inputs'][0]['sha256']='changed data'
	assert not m.s7_manifests_compatible(old,changed)
	unknown=copy.deepcopy(old);unknown['code']['abm_numeric_sha256']='unreviewed legacy source'
	assert not m.s7_manifests_compatible(unknown,new)
	assert m.s7_resume_code_hash(source.replace('a.s7_prior_strength < 0','a.s7_prior_strength < 1'))!=contract

def test_rare_continuous_signal_survives():
    x = np.zeros((10000, 1)); x[:30, 0] = 20.
    p = m.S7TailTransform().fit(x)
    z = p.transform(x)
    assert np.var(z['bulk']) == 0
    assert np.var(z['tail_up']) > 0
    assert np.all(z['tail_up'][:30] > 0) and np.all(z['tail_up'][30:] == 0)
    assert p.audit_.rare_signal_preserved_in_tail.iloc[0]

def test_rare_binary_not_winsorized():
    x = np.zeros((10000, 1)); x[:30] = 1
    p = m.S7TailTransform().fit(x); z = p.transform(x)
    np.testing.assert_array_equal(x, z['bulk'])
    assert p.feature_types_ == ('binary',)

def test_explicit_dosage_preserved():
    x = np.array([[0], [.1], [1.], [1.8], [2.]])
    p = m.S7TailTransform().fit(x, ['dosage'])
    np.testing.assert_allclose(p.transform(x)['bulk'][:, 0], x[:, 0]/2)

def test_training_constant_and_missing_dropped():
    x = np.c_[np.arange(10.), np.ones(10), np.full(10, np.nan)]
    p = m.S7TailTransform().fit(x)
    assert p.keep_.tolist() == [0]
    assert p.audit_.status.tolist() == ['retained','constant_raw','all_missing']

@pytest.mark.parametrize('x', [np.zeros((4, 1)), np.full((4, 1), np.nan)])
def test_no_usable_features_fails(x):
    with pytest.raises(ValueError): m.S7TailTransform().fit(x)

def test_projection_cannot_refit():
    rng = np.random.default_rng(5); x = rng.normal(size=(500, 3))
    p = m.S7TailTransform().fit(x); before = pickle.dumps(p)
    p.transform(np.full((5, 3), 1e10))
    assert before == pickle.dumps(p)

def test_transform_row_batch_invariance():
    x = np.random.default_rng(4).normal(size=(100, 4))
    p = m.S7TailTransform().fit(x[:60]); a = p.transform(x[60:])
    b = p.transform(x[60:][::-1])
    for name in a: np.testing.assert_array_equal(a[name], b[name][::-1])
    for name in a: np.testing.assert_array_equal(a[name][:1], p.transform(x[60:61])[name])

def test_missing_is_not_extreme():
    x = np.arange(100.)[:, None]; p = m.S7TailTransform().fit(x)
    z = p.transform(np.array([[np.nan], [np.inf], [-np.inf]]))
    assert not z['observed'].any()
    assert not z['tail_up'].any() and not z['tail_down'].any()

def test_lower_and_upper_tails_have_separate_identity():
    p = m.S7TailTransform().fit(np.arange(100.)[:, None])
    z = p.transform(np.array([[-100.], [200.]]))
    assert z['tail_down'][0, 0] > 0 and z['tail_up'][0, 0] == 0
    assert z['tail_up'][1, 0] > 0 and z['tail_down'][1, 0] == 0

def test_numerical_tail_cap_is_flagged():
    p = m.S7TailTransform().fit(np.arange(100.)[:, None])
    z = p.transform(np.array([[-1e300], [1e300]]))
    assert z['tail_cap_hit'].all()
    assert np.isfinite(z['tail_up']).all() and np.isfinite(z['tail_down']).all()

def test_linear_comparators_get_same_information():
    x = np.arange(80.).reshape(20, 4)
    p = m.S7TailTransform().fit(x)
    assert p.linear_design(x).shape == (20, 16)

def test_wrong_source_feature_count_rejected():
    p = m.S7TailTransform().fit(np.arange(10.)[:,None])
    with pytest.raises(ValueError): p.transform(np.ones((4, 2)))

def test_binary_schema_violation_fails():
    p = m.S7TailTransform().fit(np.array([[0.],[1.]]))
    with pytest.raises(ValueError): p.transform(np.array([[2.]]))

def borrow(d, **kw):
    d = np.array(d, dtype=float)
    return m.fixed_scale_borrow(d, np.ones_like(d), np.ones_like(d), np.full(len(d), .2),
                               radius=kw.pop('radius', 1.), **kw)

def test_fixed_kernel_can_focus_without_copy1():
    r = borrow([[0.] + [1.]*99], temperature=.25)
    assert r.weights[0,0] > .95
    assert r.borrow_fraction[0] < .05
    assert .2 < r.probability[0] < .3

def test_absolute_distance_reduces_support():
    near = borrow([[.1, .2, .3]])
    far = borrow([[10., 20., 30.]])
    assert far.support[0] < near.support[0]
    assert far.borrow_fraction[0] < near.borrow_fraction[0]
    assert abs(far.probability[0]-.2) < 1e-10

@pytest.mark.parametrize('mode', ['empty','masked','zero_weight','infinite'])
def test_no_donors_falls_back(mode):
    d = np.ones((2, 0 if mode == 'empty' else 3))
    w = np.zeros_like(d) if mode == 'zero_weight' else np.ones_like(d)
    if mode == 'infinite': d[:] = np.inf
    allowed = np.zeros_like(d, dtype=bool) if mode == 'masked' else None
    r = m.fixed_scale_borrow(d, np.ones_like(d), w, np.array([.1,.7]), radius=1, allowed=allowed)
    np.testing.assert_array_equal(r.probability, [.1,.7])
    assert not r.borrow_fraction.any() and not r.weights.any()

def test_related_donors_do_not_inflate_group_ess():
    r = borrow([[0.,0.,0.,0.]], donor_groups=np.array([['A','A','B','B']]))
    np.testing.assert_allclose(r.effective_donors, [4.])
    np.testing.assert_allclose(r.effective_groups, [2.])
    assert r.borrow_fraction[0] < borrow([[0.,0.,0.,0.]]).borrow_fraction[0]

def test_donor_order_invariance():
    d = np.array([[.1,.8,.3]]); y = np.array([[1.,0.,1.]]); w = np.array([[2.,1.,3.]])
    a = m.fixed_scale_borrow(d,y,w,np.array([.3]),radius=.7)
    b = m.fixed_scale_borrow(d[:,::-1],y[:,::-1],w[:,::-1],np.array([.3]),radius=.7)
    np.testing.assert_allclose(a.probability,b.probability,rtol=0,atol=1e-15)
    np.testing.assert_allclose(a.weights,b.weights[:,::-1],rtol=0,atol=1e-15)

def test_masked_label_change_has_no_effect():
    d = np.array([[.1,.2]]); y = np.array([[1.,0.]])
    kw = dict(distances=d, donor_ipcw=np.ones_like(d), prior=np.array([.2]), radius=1., allowed=np.array([[False,True]]))
    a=m.fixed_scale_borrow(labels=y,**kw); y[0,0]=0; b=m.fixed_scale_borrow(labels=y,**kw)
    np.testing.assert_array_equal(a.probability,b.probability)

def test_invalid_distance_rejected():
    with pytest.raises(ValueError): borrow([[-1.,1.]])
    with pytest.raises(ValueError): borrow([[np.nan,1.]])

def metadata():
    return pd.DataFrame({'eid':['00001','00002'], 'age':[50,60]})

def test_projection_family_file_without_internal_column():
    p=metadata(); mapping=pd.DataFrame({'eid':['00002','00001','extra'], 'group':['002','001','003']})
    z=m.attach_projection_groups(p,mapping=mapping)
    assert z.eid.tolist() == ['00001','00002']
    assert z['.le8_family'].tolist() == ['001','002']

@pytest.mark.parametrize('mode', ['missing','duplicate','conflict'])
def test_bad_family_files_fail(mode):
    p=metadata(); mp=pd.DataFrame({'eid':['00001','00002'], 'group':['A','B']})
    kw={}
    if mode=='missing': mp=mp.iloc[:1]
    elif mode=='duplicate': mp=pd.concat([mp,mp.iloc[:1]])
    else: p['family']=['A','C']; kw['source_group_col']='family'
    with pytest.raises(ValueError): m.attach_projection_groups(p,mapping=mp,**kw)

def test_projection_requires_explicit_family_information():
    with pytest.raises(ValueError): m.attach_projection_groups(metadata())
    z=m.attach_projection_groups(metadata(),allow_individual_fallback=True)
    assert z['.le8_family'].tolist() == metadata().eid.tolist()

def test_noncontext_development_person_is_blocked():
    with pytest.raises(ValueError): m.assert_external_roster(['train_not_context'],['F9'],{'train_not_context'},{'F1'})

def test_development_family_is_blocked():
    with pytest.raises(ValueError): m.assert_external_roster(['new'],['F1'],{'old'},{'F1'})

def test_external_person_and_family_pass():
    m.assert_external_roster(['new'],['new_family'],{'old'},{'old_family'})

def test_nonstring_id_rejected_instead_of_losing_zeros():
    p=metadata();p.eid=[1,2]
    with pytest.raises(ValueError): m.attach_projection_groups(p,allow_individual_fallback=True)

@pytest.mark.parametrize('device',['cpu','cuda'])
def test_real_supervised_qk_gradients_and_cpu_restore(device, tmp_path):
	if device=='cuda' and not m.torch.cuda.is_available():
		pytest.skip('Real CUDA unavailable; never substitute CPU and call it CUDA')
	a=config('--device',device,'--epochs','2','--batch-size','128','--s7-reconstruction-weight','0')
	assert a.s7_reconstruction_weight==0.;a.group_col=''
	rng=np.random.default_rng(4);n=500
	x=rng.normal(size=(n,8));p=pd.DataFrame(dict(eid=[f'device-{i}' for i in range(n)],age=rng.uniform(40,70,n),sex=rng.integers(0,2,n)))
	time=rng.exponential(8,n);event=(time<9).astype(int);time=np.minimum(time,9)
	data=m.S7FitData(x,p,time,event,m.s7_schema([f'F{i}' for i in range(8)],a),'build')
	forward_devices=[]
	original_forward=m.S7RowEncoder.forward
	def checked_forward(network,*args,**kwargs):
		forward_devices.append(next(network.parameters()).device.type)
		assert all(v.device.type==device for v in kwargs.values() if m.torch.is_tensor(v))
		return original_forward(network,*args,**kwargs)
	with threadpool_limits(2),mock.patch.object(m.S7RowEncoder,'forward',checked_forward):
		m.torch.set_num_threads(2)
		learner=m.S7AttentionLearner().fit(data.take(np.arange(400)),data.take(np.arange(400,500),'tune_model'),np.ones(400),vars(a))
		assert {'query.weight','key.weight'}<=learner.gradient_names
		assert any(name.startswith('context.') for name in learner.gradient_names)
		assert learner.fit_status['selected_tensor_changed_count']>0
		assert learner.fit_status['supervised_steps']>=2
		assert learner.fit_status['actual_device']==device
		assert learner.device==device
		assert forward_devices and set(forward_devices)=={device}
		forward_devices.clear()
		full=learner.predict(x[400:],p.iloc[400:]).probability
		assert forward_devices and set(forward_devices)=={device}
		# External prediction has finished on the requested device, then idles on CPU.
		assert next(learner.network.parameters()).device.type=='cpu'
		assert learner.device==device
	with threadpool_limits(2):
		restored=m.S7AttentionLearner.restore_state(learner.export_state(),'cpu')
		one=restored.predict(x[400:401],p.iloc[400:401]).probability
		np.testing.assert_allclose(one,full[:1],atol=1e-6,rtol=0)
		path=tmp_path/'learner.joblib';m.joblib.dump(learner,path)
		portable=m.joblib.load(path)
		assert portable.device=='cpu'
		assert all(t.device.type=='cpu' for t in portable.network.state_dict().values())
		assert m.s7_state_hash(portable.network)==m.s7_state_hash(learner.network)
		frozen=m.s7_state_hash(learner)
		learner.device='cpu'
		assert m.s7_state_hash(learner)==frozen
		learner.device=device
		np.testing.assert_allclose(portable.predict(x[400:401],p.iloc[400:401]).probability,full[:1],atol=1e-6,rtol=0)
	print(dict(device=device,**learner.fit_status))

def test_default_device_refuses_cpu_fallback_before_inputs():
	a=m.configure(m.reference_parser().parse_args(['--demo','--preflight']))
	assert a.abm_design=='selective_attention' and a.device=='cuda'
	with mock.patch.object(m.torch.cuda,'is_available',return_value=False),mock.patch.object(m,'preflight_inputs') as inputs:
		with pytest.raises(RuntimeError,match='refusing CPU fallback'):
			m.s7_main(a)
		inputs.assert_not_called()

def test_cuda_kernel_failure_refuses_cpu_fallback():
	with mock.patch.object(m.torch.cuda,'is_available',return_value=True),mock.patch.object(m.torch,'ones',side_effect=RuntimeError('kernel unavailable')):
		with pytest.raises(RuntimeError,match='CUDA forward/backward failed.*refusing CPU fallback'):
			m.check_execution_device('cuda')

@pytest.mark.parametrize('scenario', ['S01', 'S02', 'S03', 'S04', 'S05'])
def test_structural_simulation(scenario):
	root=Path(os.getenv('LE8_S7_SIMULATION_DIR') or tempfile.mkdtemp(prefix='s7-sim-',dir='/tmp'))/scenario
	root.mkdir(parents=True,exist_ok=True)
	rng=np.random.default_rng(2026+int(scenario[1:]))
	n=10000 if scenario=='S01' else 3000
	x=rng.normal(size=(n,8));z=rng.integers(0,2,n)
	p=pd.DataFrame(dict(eid=[f'{scenario}-{i:05d}' for i in range(n)],family=[f'{scenario}-family-{i//2}' for i in range(n)],
		age=rng.uniform(40,70,n),sex=rng.integers(0,2,n)))
	a=config('--epochs','6','--patience','3','--width','16','--tokens','4','--batch-size','128',
		'--tree-estimators','40','--s7-radius-ref-n','100','--c-grid','.01,.1,1','--s7-bandwidth-grid','1')
	if scenario=='S01':
		rare=rng.random(n)<.003;x[:,5]=20*rare
		eta=-3+3*rare+.3*x[:,0];oracle_features=np.c_[x,rare.astype(float)]
	elif scenario in ('S02','S03','S04'):
		x[:,:4]=rng.integers(0,2,(n,4))
		interaction=(x[:,0]==x[:,1])
		if scenario!='S02':interaction=np.where(z==1,interaction,x[:,2]==x[:,3])
		eta=m.logit(np.where(interaction,.6,.05))
		oracle_features=np.c_[x,interaction.astype(float)]
		if scenario=='S03':x=np.c_[x,z]
	else:
		plate=rng.random(n)<np.where(np.arange(n)<int(.8*n),.003,.3)
		x[:,5]=20*plate;p['plate']=plate.astype(str)
		eta=-1+.5*x[:,0];oracle_features=x[:,[0]]
		a.residualize='plate';a.categorical='sex,plate'
	prob=m.expit(eta);failure=rng.exponential(5/-np.log1p(-prob));censor=rng.uniform(6,9,n)
	time=np.minimum(failure,censor);event=(failure<=censor).astype(int)
	train=np.arange(int(.6*n));valid=np.arange(int(.6*n),int(.8*n));test=np.arange(int(.8*n),n)
	names=[f'F{i}' for i in range(x.shape[1])]
	schema=m.s7_schema(names,a)
	data=m.S7FitData(x,m.s7_metadata(p,a),time,event,schema,'build')
	# Oracle-only information is saved separately and never passed to actual learners.
	pd.DataFrame(x,columns=names).assign(eid=p.eid).to_csv(root/'model_inputs.csv',index=False)
	pd.DataFrame(dict(eid=p.eid,oracle_probability=prob,latent_Z=z,time=time,event=event)).to_csv(root/'oracle_truth.csv',index=False)
	predictions={'oracle_truth':prob[test]};status=[]
	with threadpool_limits(2):
		m.torch.set_num_threads(2)
		for name,learner in [('elasticnet',m.S7SklearnLearner()),('hist',m.S7SklearnLearner('hist')),
			('attention_direct',m.S7AttentionLearner('direct')),('attention_retrieval',m.S7AttentionLearner('retrieval'))]:
			try:
				learner.fit(data.take(train),data.take(valid,'tune_model'),np.ones(len(train)),vars(a))
				predictions[name]=learner.predict(x[test],data.metadata.iloc[test]).probability
				status.append(dict(model=name,status='completed',**learner.fit_status))
				if scenario=='S01':
					ix=learner.prep.feature_names.index('F5');batch=learner.prep.transform(x[train],data.metadata.iloc[train])
					assert np.any(batch.tail_up[:,ix]>0)
				if scenario=='S05':assert learner.prep.technical is not None
				learner.prep.audit.assign(model=name).to_csv(root/(name+'_preprocessing.csv'),index=False)
			except Exception as exc:
				status.append(dict(model=name,status='failed',error=str(exc)))
				pd.DataFrame(status).to_csv(root/'fit_status.csv',index=False)
				raise
		oracle=m.S7SklearnLearner()
		oracle_data=m.S7FitData(oracle_features,data.metadata,time,event,m.s7_schema([f'O{i}' for i in range(oracle_features.shape[1])],a),'build')
		oracle.fit(oracle_data.take(train),oracle_data.take(valid,'tune_model'),np.ones(len(train)),{**vars(a),'c_grid':[.1,1,10]})
		predictions['oracle_regression']=oracle.predict(oracle_features[test],data.metadata.iloc[test]).probability
		status.append(dict(model='oracle_regression',status='completed',**oracle.fit_status))
		km=m.S6Censoring().fit(data.take(train).survival(),5)
		groups=p.family.to_numpy()[test]
		metrics=pd.DataFrame([dict(scenario=scenario,model=model,**m.s7_metric(time[test],event[test],values,groups,km)) for model,values in predictions.items()])
		pairs=[(name,'oracle_truth') for name in predictions if name!='oracle_truth']
		pairs += [('attention_retrieval','elasticnet'),('attention_retrieval','attention_direct')]
		ci=m.s7_paired_bootstrap(time[test],event[test],predictions,np.ones(len(test),bool),groups,km,pairs,30,a.seed)
	pd.DataFrame(status).to_csv(root/'fit_status.csv',index=False);metrics.to_csv(root/'metrics.csv',index=False);ci.to_csv(root/'paired_intervals.csv',index=False)
	(root/'interpretation.json').write_text(json.dumps(dict(scenario=scenario,synthetic=True,N=n,clinical_claim=False,
		neural_epochs=6,bootstrap=30,uncertainty='conditional_on_frozen_models; exploratory Monte Carlo interval',
		rare_count=int((x[:,5]>0).sum()) if scenario=='S01' else None,
		latent_Z_available_to_model=scenario=='S03',subgroup_mechanisms_identified=False),indent=2))
	assert np.isfinite(metrics.loc[metrics.model.isin(['elasticnet','attention_direct','attention_retrieval']),'Brier_IPCW']).all()
	# This is an implementation/capability experiment, never an assertion that neural must win.
	print(metrics[['scenario','model','AUC_IPCW','Uno_C_horizon','Brier_IPCW']].to_string(index=False))

def test_two_level_074_is_not_individual_predictability():
	event=np.r_[np.ones(74),np.zeros(26),np.ones(26),np.zeros(74)].astype(int)
	time=np.where(event,1.,8.);risk=np.repeat([.74,.26],100)
	km=m.S6Censoring().fit(pd.DataFrame(dict(time=time,event=event)),5)
	assert m.s7_uno(time,event,risk,km)['Uno_C_horizon']==pytest.approx(.74)
	for part in (slice(0,100),slice(100,200)):
		assert m.s7_uno(time[part],event[part],risk[part],km)['Uno_C_horizon']==pytest.approx(.5)

def test_public_cpu_frozen_roundtrip_and_final():
	root=Path(os.getenv('LE8_S7_PUBLIC_DIR') or tempfile.mkdtemp(prefix='s7-public-',dir='/tmp'))
	inputs=Path(os.getenv('LE8_S7_INPUT_DIR') or root.parent/(root.name+'-synthetic-input'))
	root.mkdir(parents=True,exist_ok=True)
	if not (inputs/'phe.csv').exists():
		subprocess.run([sys.executable,str(ROOT/'validation/c1_abm.py'),'--smoke-data','--out',str(inputs)],check=True)
	out=root/'cvd_cad/prot/c1_correlate/abm_selective_attention'
	args=['--abm-design','selective_attention','--phe-file',str(inputs/'phe.csv'),'--omics-file',str(inputs/'prot.csv'),
		'--covariates','age,sex','--categorical','sex','--residualize','','--horizon','5','--tree','hist','--tree-estimators','10',
		'--pretrain-epochs','1','--epochs','2','--patience','2','--tokens','8','--width','16','--heads','2','--layers','1',
		'--batch-size','64','--device','cpu','--cores','2','--selective-folds','3','--selective-min-events','5',
		'--s7-neighbors-per-view','16','--bootstrap','30']
	base=[str(ROOT/'le8.sh'),'c1_abm','--Y','cvd_cad','--biom','prot','--analysis-root',str(root),'--abm-backend','reference']
	train=base+['--group-file',str(inputs/'groups.csv'),'--outer-roster',str(inputs/'outer.csv'),
		'--end-date','2020-01-01','--seed','2026','--abm-args',shlex.join(args)]
	def run(command,name):
		result=subprocess.run(command,cwd=ROOT,env={**os.environ,'ABM_PYTHON':sys.executable,'PYTHONDONTWRITEBYTECODE':'1'},text=True,capture_output=True)
		(root.parent/(root.name+'-'+name+'.log')).write_text(result.stdout+result.stderr)
		assert result.returncode==0,result.stdout[-4000:]+result.stderr[-4000:]
		for path in re.findall(r'Temporary tables and execution files: (\S+)',result.stdout):
			assert not Path(path).exists(),'Successful temporary workspace should be removed'
		return result
	run(train,'train')
	assert (out/'model_bundle.joblib').is_file()
	sha=lambda:hashlib.sha256((out/'model_bundle.joblib').read_bytes()).hexdigest()
	before=sha();cached=run(train,'cache')
	assert 'SKIP' in cached.stdout and sha()==before
	assert not list(out.glob('*.csv'))
	assert not (out/'test_individuals.rds').exists()
	from zipfile import ZipFile
	assert any('test_individuals.csv' in json.loads(ZipFile(p).read('le8/manifest.json')).get('files', {}) for p in out.glob('*.xlsx') if 'le8/manifest.json' in ZipFile(p).namelist())
	external=root.parent/(root.name+'-external');external.mkdir(exist_ok=True)
	p=pd.read_csv(inputs/'phe.csv').iloc[:12][['eid','age','sex']].copy()
	x=pd.read_csv(inputs/'prot.csv').iloc[:12].copy()
	p.eid=[f'external-{i:04d}' for i in range(len(p))];x.eid=p.eid
	p.to_csv(external/'phe.csv',index=False);x.to_csv(external/'omics.csv',index=False)
	pd.DataFrame(dict(eid=p.eid,group=[f'external-family-{i//2}' for i in range(len(p))])).to_csv(external/'groups.csv',index=False)
	project=['project','--run-dir',str(out),'--phe-file',str(external/'phe.csv'),'--omics-file',str(external/'omics.csv'),
		'--output',str(external/'prediction.csv'),'--projection-group-file',str(external/'groups.csv'),
		'--projection-group-namespace','ukb','--device','cpu','--cores','2']
	# Remove access to the synthetic training paths; only the frozen artifact may be needed.
	hidden=inputs.with_name(inputs.name+'-temporarily-unavailable');assert not hidden.exists()
	inputs.rename(hidden)
	try:
		run(base+['--abm-args',shlex.join(project)],'project')
		one=pd.read_csv(external/'prediction.csv').set_index('eid')
		p.iloc[::-1].to_csv(external/'phe.csv',index=False);x.iloc[::-1].to_csv(external/'omics.csv',index=False)
		run(base+['--abm-args',shlex.join(project)],'project-reordered')
		two=pd.read_csv(external/'prediction.csv').set_index('eid').loc[one.index]
		np.testing.assert_allclose(one.policy_risk,two.policy_risk,atol=1e-6,rtol=0)
		assert (one.validation_status=='external_validation').all()
		run(base+['--abm-args',shlex.join(['evaluate','--run-dir',str(out),'--device','cpu','--cores','2'])],'evaluate')
	finally: hidden.rename(inputs)
	assert sha()==before
	# A real old S6 fit supplies the side-by-side legacy source for Final.

	from threadpoolctl import threadpool_limits
	legacy=root/'cvd_cad/prot/c1_correlate/abm_reference'
	if not (legacy/'model_bundle.joblib').exists():
		old_a=config();old_a.abm_design='selective';old_a.trait='cvd_cad'
		xold,pold,names=synthetic();pold['split']=m.s6_split(pold,old_a)
		with threadpool_limits(2):
			m.s6_fit(xold,pold,names,old_a,legacy);m.s6_evaluate(legacy)
		subprocess.run(['Rscript',str(ROOT/'f/0.common.R'),'--tables-pack',str(legacy)],check=True,
			env={**os.environ,'LE8_TABLE_WORKSPACE':'1'},capture_output=True,text=True)
	run([str(ROOT/'le8.sh'),'final','--Y','cvd_cad','--biom','prot','--analysis-root',str(root)],'final')
	assert (root/'final/Fig9.question_ABM_attention.png').is_file()
	assert (root/'final/Fig9.question_ABM_attention.xlsx').is_file()
	# Restore a COPY to inspect exact aggregate rows, without unpacking the published root.
	inspect=root.parent/(root.name+'-restored');inspect.mkdir(exist_ok=True)
	shutil.copytree(root/'final',inspect/'final',dirs_exist_ok=True)
	shutil.copytree(root/'cvd_cad',inspect/'cvd_cad',dirs_exist_ok=True)
	subprocess.run(['Rscript',str(ROOT/'f/0.common.R'),'--tables-restore',str(inspect)],check=True,
		env={**os.environ,'LE8_TABLE_WORKSPACE':'1','LE8_TABLE_ABM_PRIVATE':'TRUE'},capture_output=True,text=True)
	metrics=pd.read_csv(inspect/'final/final.questions.abm_metrics.csv')
	registry=pd.read_csv(inspect/'final/final.questions.abm_registry.csv')
	assert set(metrics.backend)=={'reference','selective_attention'}
	assert set(registry.model_id)<=set(metrics.loc[metrics.backend.eq('selective_attention'),'model_id'])
	assert registry.architecture.str.contains('S7RowEncoder').any()
	assert metrics.model_id.nunique()>=10
	roles=pd.read_csv(inspect/'cvd_cad/prot/c1_correlate/abm_selective_attention/development_roster.csv')
	assert roles.groupby('.le8_family').role.nunique().max()==1
	assert 'tune_gate' in set(roles.role)
	(root.parent/(root.name+'-verification.json')).write_text(json.dumps(dict(public_cpu=True,cache=True,
		individual_XLSX_roundtrip=True,outcome_free_external=True,training_inputs_unavailable=True,model_sha256=sha(),
		final_registered_models=int(registry.model_id.nunique()),synthetic_only=True),indent=2))

def legacy_config(**changes):
    p=argparse.ArgumentParser();m.s6_options(p);d=vars(p.parse_args([]))
    d.update(seed=2026,trait="synthetic",biom="prot",id_col="eid",group_col="family",covariates="age,sex",
      residualize="",categorical="sex",transform="none",feature_missing=.2,sample_missing=.2,
      horizon=5,min_censor_survival=.05,teacher_c=.1,max_iter=2000,tree="hist",tree_estimators=20,
      cores=2,c_grid=[.01,.1],bootstrap=20,split_file="",module_file="",resume=False,
      selective_folds=3,selective_gate_trees=20,selective_components=6,selective_neighbors=30,
      selective_min_events=5,selective_coverages="0.4,0.6,1.0",shuffle_development_outcomes=False)
    d.update(changes)
    return m.s6_validate(types.SimpleNamespace(**d))

def legacy_data(n=2400,seed=315):
    rng=np.random.default_rng(seed)
    x=rng.normal(size=(n,12)).astype("float32")
    age=rng.uniform(40,70,n);sex=rng.integers(0,2,n)
    rate=.025*np.exp(.035*(age-55)+.5*sex+.9*x[:,0]-.6*x[:,1])
    t=rng.exponential(1/rate);c=rng.uniform(6,10,n)
    p=pd.DataFrame(dict(eid=["p"+str(i) for i in range(n)],family=["f"+str(i//2) for i in range(n)],
         age=age,sex=sex,time=np.minimum(t,c),event=(t<=c).astype(int)))
    return x,p,["X"+str(i) for i in range(x.shape[1])]

@threadpool_limits.wrap(limits=2)
def test_legacy_selective_training():
    passed=[]
    def done(name):passed.append(name);print("PASS "+name,flush=True)
    a=legacy_config();x,p,features=legacy_data()
    groups=m.s6_groups(p,a)
    fold=m.s6_group_folds(groups,3,90)
    for g in np.unique(groups):assert len(set(fold[groups==g]))==1
    shuffled=np.random.default_rng(2).permutation(len(p))
    assert np.array_equal(fold[shuffled],m.s6_group_folds(groups[shuffled],3,90))
    done("outcome-blind, family-preserving and row-order-stable folds")
    try:legacy_config(covariates="age,event")
    except ValueError:pass
    else:raise AssertionError("Outcome predictor was accepted")
    done("outcome/follow-up covariates rejected")
    ids=np.array([f"a{i}" for i in range(100)])
    score=np.zeros(100);support=np.linspace(0,1,100)
    rule=m.S6CoverageRule().fit(score,support,ids,.6,.99,7)
    assert rule.apply(score,support,ids).sum()==60
    full=m.S6CoverageRule().fit(score,support,ids,1.,.1,7)
    assert full.apply(score,support,ids).all()
    repeat=rule.apply(score[shuffled[:100]%100],support[shuffled[:100]%100],ids[shuffled[:100]%100])
    assert np.array_equal(repeat,rule.apply(score,support,ids)[shuffled[:100]%100])
    done("constant-score ties and true 100-percent coverage endpoint")
    # Pure probability borrowing must fall back instead of using a family donor.
    small_a=legacy_config(selective_neighbors=5)
    bank=m.S6Neighborhood().fit(x[:100],np.tile([0,1],50),np.ones(100),np.array(["same"]*100),small_a,3,features)
    prob,diag=bank.predict(x[100:105],np.full(5,.2),np.array(["same"]*5))
    assert np.allclose(prob,.2) and np.all(diag[:,2]==0)
    done("same-family donors excluded; empty neighborhood uses prior")
    prep=m.S6Preprocessor(.2,[],["sex"],"none").fit(x[:1000],p.iloc[:1000])
    mean=prep.mean.copy();median=prep.median.copy()
    prep.transform(x[1000:1010]*1e6,p.iloc[1000:1010])
    assert np.array_equal(mean,prep.mean) and np.array_equal(median,prep.median)
    done("test transformations cannot change fitted preprocessing")
    # Nontrivial censoring and sparse plate correction exercise UKB-specific paths.
    mixed=p.copy();rr=np.random.default_rng(55)
    early=rr.choice(len(p),len(p)//4,replace=False)
    mixed.loc[early,"time"]=rr.uniform(.1,3.,len(early));mixed.loc[early,"event"]=0
    censor=m.S6Censoring().fit(mixed,a.horizon,a.min_censor_survival)
    yy,ww=censor.labels_weights(mixed)
    assert (ww[early]==0).all() and (yy[early]==0).all()
    assert ww.max()>1 and np.isfinite(ww).all()
    assert censor.at([a.horizon])[0]<1
    done("pre-horizon censoring is not mislabeled as an observed control")
    technical=p.iloc[:600].copy();technical["plate"]=np.arange(600)%60
    z=x[:600].copy();z[np.arange(60),np.arange(60)%12]=np.nan
    corrected=m.S6Preprocessor(.2,["plate"],["sex","plate"],"none").fit(z,technical)
    zz,observed=corrected.transform(z,technical)
    assert np.isfinite(zz).all() and (~observed).sum()==60
    assert m.sparse.issparse(corrected.design.transform(technical))
    done("sparse technical residualization and missing molecular measurements")
    # Uno-style implementation checked against its explicit comparable-pair formula.
    km=m.S6Censoring().fit(p,5,.05)
    rng=np.random.default_rng(22);t=rng.integers(1,9,45).astype(float);e=rng.integers(0,2,45);risk=rng.integers(0,5,45)/5
    num=den=0.
    for i in range(len(t)):
        if e[i] and t[i]<=km.horizon:
            wi=1/km.at([t[i]],left=True)[0]**2
            for j in range(len(t)):
                if t[j]>t[i]:
                    den+=wi;num+=wi*(float(risk[i]>risk[j])+.5*float(risk[i]==risk[j]))
    assert abs(m.s6_uno(t,e,risk,km)-num/den)<1e-12
    assert abs(m.s6_uno(t,e,np.ones(len(t)),km)-.5)<1e-12
    done("IPCW concordance equals naive pair calculation; ties score 0.5")
    with tempfile.TemporaryDirectory(prefix="c1-selective-check-") as tmp:
        tmp=Path(tmp)
        module_file=tmp/"modules.tsv"
        pd.DataFrame({"feature":features[:6],"module":["A"]*3+["B"]*3}).to_csv(module_file,sep="\t",index=False)
        module_a=legacy_config(module_file=str(module_file))
        bank=m.S6Neighborhood().fit(x[:100],np.tile([0,1],50),np.ones(100),groups[:100],module_a,4,features)
        assert bank.transform(x[:3]).shape[1]==8 and bank.module_names==["A","B"]
        done("predeclared molecular modules included in geometry")
        # Flip the outcomes of one entire recipient fold. Its own raw gate scores,
        # hard masks and soft weights must be bit-for-bit unchanged.
        cf=legacy_config(tree="none",selective_gate_trees=15)
        base_x,base_p=x[:1800],p.iloc[:1800].copy()
        held=m.s6_group_folds(m.s6_groups(base_p,cf),cf.selective_folds,cf.seed+1200)==0
        _,v1=m.s6_crossfit_gate(base_x,base_p,features,cf)
        changed=base_p.copy();changed.loc[held,"event"]=1-changed.loc[held,"event"]
        _,v2=m.s6_crossfit_gate(base_x,changed,features,cf)
        for key in ["omics_gain","local_gain","training_selected","training_soft"]:
            assert np.array_equal(v1[key][held],v2[key][held]),key
        done("own-label poisoning cannot change honest fold-local training weights")
        keep,soft=m.s6_train_weights(v1,base_p.eid.to_numpy(),cf)
        rnd,rw=m.s6_matched_random(keep,soft,v1["OOF_clinical"],base_p.eid.to_numpy(),2)
        assert rnd.sum()==keep.sum() and np.allclose(np.sort(rw),np.sort(soft))
        assert np.min(soft)>=cf.selective_weight_floor and np.max(soft)<=1
        done("matched random sample size/weight distribution; bounded soft weights")
        # End-to-end: all fitted arms, calibration, frozen prediction, metrics, PNGs.
        p=p.copy();p["split"]=m.s6_split(p,a)
        out=tmp/"first";b=m.s6_fit(x,p,features,a,out);result=m.s6_evaluate(out)
        table=pd.read_csv(out/"test_individuals.csv")
        for _,sub in result.groupby("subset"):
            assert sub.n.nunique()==1 and sub.cases.nunique()==1
        assert {"AUC_IPCW","Uno_C_horizon","Brier_IPCW"}<=set(result)
        assert (out/"Fig_selective_training.png").is_file()
        done("end-to-end survival fitting, calibration, same-mask comparisons and figures")
        # Real inference does not need outcome/follow-up columns.
        test=np.flatnonzero(p.split.eq("test"))
        projection=m.s6_predict(b,x[test],p.iloc[test].drop(columns=["time","event"]))
        assert np.allclose(projection[b["primary"]],table[b["primary"]],atol=1e-12)
        done("serialized model reproduces outcome-free inference")
        # Refit with a permuted TEST outcome file and unchanged frozen split.
        p2=p.copy();perm=np.random.default_rng(200).permutation(test)
        p2.loc[test,["time","event"]]=p2.loc[perm,["time","event"]].to_numpy()
        b2=m.s6_fit(x,p2,features,a,tmp/"poisoned_test")
        t2=pd.read_csv(tmp/"poisoned_test"/"test_individuals.csv")
        for col in list(b["model"].calibrators)+["omics_gain","local_gain","supported_match","prediction_released"]:
            assert np.allclose(table[col],t2[col],equal_nan=True,atol=1e-12),col
        done("test-label permutation leaves every prediction and selection mask unchanged")
        before=hashlib.sha256((out/"model_bundle.joblib").read_bytes()).hexdigest()
        def forbidden_fit(*_a,**_k):raise AssertionError("Evaluation attempted model fitting")
        with patch.object(m.S6Calibrator,"fit",forbidden_fit),patch.object(m.S6Preprocessor,"fit",forbidden_fit),patch.object(m.S6GainGate,"fit",forbidden_fit):
            m.s6_evaluate(out)
        assert before==hashlib.sha256((out/"model_bundle.joblib").read_bytes()).hexdigest()
        done("evaluation is fit-free and cannot mutate frozen model")
        if importlib.util.find_spec("lightgbm"):
            la=legacy_config(tree="lightgbm")
            learner=m.s6_fit_tree(m.s6_make_tree(la,4),x[:300],(x[:300,0]>0).astype(int),np.ones(300))
            assert np.all(np.isfinite(learner.predict_proba(x[300:310])))
            done("installed LightGBM branch trains and predicts")
        else:print("SKIP LightGBM not installed")
    assert m.reference_parser().parse_args([]).abm_design=="selective_attention"
    done("actual integrated host import and selective CLI")
    print(f"\n{len(passed)} checks passed. Synthetic data only; no UKB training and no R model execution.")

def make_smoke_data(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n',type=int,default=6000)
    ap.add_argument('--features',type=int,default=24)
    ap.add_argument('--seed',type=int,default=2026)
    a=ap.parse_args(argv)
    if a.n < 1000 or a.n%2 or a.features < 8:
        ap.error('Use an even N >= 1000 and >= 8 features')
    a.out.mkdir(parents=True,exist_ok=True)
    names=['phe.csv','prot.csv','groups.csv','outer.csv','SIMULATION_NOTICE.json']
    if any((a.out/n).exists() for n in names):
        raise FileExistsError('Refusing to overwrite prior inputs; choose a new directory')
    rng=np.random.default_rng(a.seed)
    ids=[f'p{i:07d}' for i in range(a.n)]
    group=np.array([f'fam{i//2:06d}' for i in range(a.n)])
    age=rng.uniform(40,70,a.n); sex=rng.integers(0,2,a.n)
    x=rng.normal(size=(a.n,a.features))
    rare=rng.random(a.n)<.003
    x[:,5]=20.*rare
    x[:,6]=rng.integers(0,2,a.n);x[:,7]=rng.integers(0,2,a.n)
    family_effect=np.repeat(rng.normal(0,.15,a.n//2),2)
    # Enough events for role/interface smoke; not a realistic CAD prevalence model.
    rate=.10*np.exp(.03*(age-55)+.25*sex+.35*x[:,0]-.25*x[:,1]+.2*(x[:,6]==x[:,7])+family_effect+.5*rare)
    failure=rng.exponential(1/rate);censor=rng.uniform(6,9,a.n)
    event=failure<=censor
    baseline=pd.Timestamp('2010-01-01')
    end=baseline+pd.to_timedelta(np.rint(censor*365.25).astype(int),unit='D')
    # Only convert observed diagnoses: rare enormous hypothetical times are not material.
    diagnosis=np.full(a.n,'',dtype=object)
    days=np.maximum(1,np.rint(failure[event]*365.25).astype(int))
    diagnosis[event]=(baseline+pd.to_timedelta(days,unit='D')).strftime('%Y-%m-%d')
    phe=pd.DataFrame(dict(eid=ids,age=age,sex=sex,date_attend='2010-01-01',
        date_lost=end.strftime('%Y-%m-%d'),date_death='',fod_icd10_cvd_cad=diagnosis))
    omics=pd.DataFrame(x,columns=[f'F{i:03d}' for i in range(a.features)])
    omics.insert(0,'eid',ids)
    # Fold definition depends only on family IDs and seed, not failure outcomes.
    levels=sorted(set(group),key=lambda v:hashlib.blake2b(f'{a.seed}:{v}'.encode(),digest_size=8).digest())
    testing=set(levels[:round(.2*len(levels))])
    roles=np.array(['test' if v in testing else 'training' for v in group])
    phe.to_csv(a.out/'phe.csv',index=False)
    omics.to_csv(a.out/'prot.csv',index=False)
    pd.DataFrame(dict(eid=ids,group=group)).to_csv(a.out/'groups.csv',index=False)
    pd.DataFrame(dict(eid=ids,role=roles)).to_csv(a.out/'outer.csv',index=False)
    notice=dict(synthetic=True,clinical_use=False,N=a.n,features=a.features,seed=a.seed,
        events_observed=int(event.sum()),events_by_5years=int(((failure<=5)&event).sum()),
        test_N=int((roles=='test').sum()),rare_count=int(rare.sum()),
        scope='Input/interface smoke only; no UKB data and no model performance claim')
    (a.out/'SIMULATION_NOTICE.json').write_text(json.dumps(notice,indent=2)+'\n')
    print(json.dumps(notice,indent=2))

class C1InputRegression(TestCase):
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

	def test_r_python_outcome_clock_parity(self):
		a=self.config('--end-date','2020-12-31')
		p=pd.DataFrame(dict(eid=[f'p{i}' for i in range(8)],birth_date='1960-01-01',date_attend=['2020-01-01']*7+['2020-12-31'],
			date_death=[None,None,None,None,'2020-06-01',None,None,None],date_lost=[None]*5+['2020-04-01',None,None],
			fod_icd10_cvd_cad=['2019-12-31','2020-01-01','2020-12-31','2021-01-01','2020-07-01','2020-04-01',None,None]))
		file=self.path/'outcomes.csv';p.to_csv(file,index=False)
		r=subprocess.run(['Rscript',str(ROOT/'validation/native_models.R'),'--outcomes',str(file)],env=dict(os.environ,DATE_FOLLOW_END=a.end_date),text=True,capture_output=True,check=True)
		data=json.loads(r.stdout);out,audit=abm.outcomes(p,a)
		self.assertEqual(data['prevalent'],out.prevalent.astype(int).tolist())
		self.assertEqual(data['event'],[int(x) if valid else None for x,valid in zip(out.event,out.endpoint_valid)])
		np.testing.assert_allclose([np.nan if x is None else x for x in data['time']],out.time.where(out.endpoint_valid).to_numpy(),equal_nan=True)


@dataclasses.dataclass
class Batch:
    bulk: np.ndarray
    tail_up: np.ndarray
    tail_down: np.ndarray
    observed: np.ndarray
    clinical: np.ndarray
    ids: np.ndarray
    groups: np.ndarray
    feature_names: list

    def design(self):
        return np.c_[self.bulk, self.tail_up, self.tail_down, (~self.observed).astype(float), self.clinical]

    def take(self, ix):
        return Batch(**{k: v if k == 'feature_names' else v[ix] for k, v in vars(self).items()})

def bank(n=1000, kind='binary'):
    x = np.zeros((n, 2))
    up = np.zeros_like(x)
    if kind == 'binary':
        x[:30, 0] = 1
    else:
        up[:30, 0] = 2
    return Batch(x, up, np.zeros_like(x), np.ones_like(x, dtype=bool), np.empty((n, 0)), np.array(['id' + str(i) for i in range(n)]), np.array(['g' + str(i) for i in range(n)]), ['a', 'b'])

def cfg():
    return dict(seed=19, s7_views=['learned', 'tail'], s7_radius_ref_n=50, s7_neighbors_per_view=8)

@pytest.mark.parametrize('kind', ['binary', 'continuous'])
def test_exact_rare_survives_zero_radius(kind):
    b = bank(kind=kind)
    z = np.zeros((len(b.ids), 2))
    g = m.S7RetrievalGeometry().fit(b, z, cfg())
    q = b.take([0])
    q.ids[:] = 'external'
    q.groups[:] = 'external'
    ix, r, ok, d = g.candidates(q, np.zeros((1, 2)))
    assert g.status['tail'] == 'exact_match_only'
    assert ok.sum() == 8 and np.all(ix[ok] < 30) and np.all(r[ok] == 0)
    assert d.tail_available.iloc[0] and (not d.learned_available.iloc[0])

@pytest.mark.parametrize('kind', ['binary', 'continuous'])
def test_exact_does_not_match_missing_or_other_assay(kind):
    b = bank(kind=kind)
    g = m.S7RetrievalGeometry().fit(b, np.zeros((len(b.ids), 2)), cfg())
    q = b.take([0])
    q.ids[:] = 'external'
    q.groups[:] = 'external'
    q.observed[:, 0] = False
    assert not g.candidates(q, np.zeros((1, 2)))[2].any()
    q.observed[:] = True
    q.bulk[:] = 0
    q.tail_up[:] = 0
    q.tail_up[:, 1] = 2
    assert not g.candidates(q, np.zeros((1, 2)))[2].any()

@pytest.mark.parametrize('exclude', ['family', 'fold'])
def test_exact_excludes_family_and_fold(exclude):
    b = bank()
    g = m.S7RetrievalGeometry().fit(b, np.zeros((len(b.ids), 2)), cfg())
    q = b.take([0])
    q.ids[:] = 'outside'
    if exclude == 'family':
        g.groups[:30] = 'F'
        q.groups[:] = 'F'
        args = {}
    else:
        q.groups[:] = 'new'
        bf = np.zeros(len(b.ids), int)
        bf[:30] = 1
        args = dict(query_folds=np.array([1]), bank_folds=bf)
    assert not g.candidates(q, np.zeros((1, 2)), **args)[2].any()

def test_constant_global_not_mislabeled_exact():
    b = bank()
    b.bulk[:] = 0
    b.tail_up[:] = 0
    g = m.S7RetrievalGeometry().fit(b, np.zeros((len(b.ids), 2)), cfg())
    assert g.status['learned'] == 'constant_geometry' and g.status['tail'] == 'unavailable'

def test_static_refresh_preserves_indexes_and_radii():
    b = bank()
    rng = np.random.default_rng(1)
    z = rng.normal(size=(len(b.ids), 2))
    g = m.S7RetrievalGeometry().fit(b, z, cfg())
    tail = g.tail_index
    radius = g.radii['tail']
    g.refresh_learned(b, z * 2)
    assert g.tail_index is tail and g.radii['tail'] == radius and (g.static_build_count == 1)
    with pytest.raises(ValueError):
        g.refresh_learned(b.take(np.arange(len(b.ids))[::-1]), z)

def test_geometry_query_batch_and_order_invariance():
    b = bank()
    z = np.random.default_rng(3).normal(size=(len(b.ids), 3))
    geo = m.S7RetrievalGeometry().fit(b, z, cfg())
    q = b.take([0, 4])
    q.ids[:] = ['x', 'y']
    q.groups[:] = ['xg', 'yg']
    one = geo.candidates(q, z[[0, 4]])
    two = geo.candidates(q.take([1, 0]), z[[4, 0]])
    np.testing.assert_array_equal(one[0], two[0][::-1])
    np.testing.assert_allclose(one[1], two[1][::-1])

def test_policy_technical_abstention_and_regular_fallback(fitted):
    b=fitted['b'];p=fitted['p'].iloc[:5].copy();x=fitted['x'][:5].copy();x[0,:]=np.nan
    table=m.s7_predict(b,x,p)
    assert not table.policy_prediction_available.iloc[0] and np.isnan(table.policy_risk.iloc[0])
    assert table.policy_source.iloc[0]=='abstain_technical_QC'
    ordinary=(~table.released)&table.policy_prediction_available
    assert ordinary.any()
    np.testing.assert_allclose(table.loc[ordinary,'policy_risk'],table.loc[ordinary,b['fallback']])
    metric=m.s7_metric(p.time.to_numpy(),p.event.to_numpy(),table.policy_risk.to_numpy(),p.family.to_numpy(),b['censoring'])
    assert metric['eligible_N']==5 and metric['abstained_N']>=1 and metric['policy_available_N']<5

def test_calibration_uses_frozen_cuts_and_preserves_artifact(fitted):
    b=fitted['b'];out=fitted['out'];before=hashlib.sha256((out/'model_bundle.joblib').read_bytes()).hexdigest()
    m.s7_evaluate(out)
    assert before==hashlib.sha256((out/'model_bundle.joblib').read_bytes()).hexdigest()
    d=pd.read_csv(out/'test_calibration.csv');g=pd.read_csv(out/'test_calibration_groups.csv')
    assert set(d.cutpoint_source)=={'calibration_fit'} and (d.policy_available_N+d.abstained_N==d.eligible_N).all()
    assert set(['all','research','released'])<=set(d.subset)
    for name,rows in g.groupby('model_id'):
        expected=set(np.asarray(b['calibration_cutpoints'][name]).tolist())
        for edge in rows.upper[np.isfinite(rows.upper)]:assert np.min(np.abs(np.asarray(list(expected))-edge))<1e-14

def test_frozen_cpu_cuda_probabilities_gate_release(fitted):
    import copy
    b=copy.deepcopy(fitted['b']);x=fitted['x'][:32];p=fitted['p'].iloc[:32]
    cpu=m.s7_predict(b,x,p)
    for learner in b['models'].values():
        if isinstance(learner,m.S7AttentionLearner):
            learner.device='cuda';learner.config['s7_retrieval_device']='cuda'
    gpu=m.s7_predict(b,x,p)
    for name in ('policy_risk',b['primary'],'gain_score'):
        np.testing.assert_allclose(cpu[name],gpu[name],atol=2e-6,rtol=0,equal_nan=True)
    for name in ('research_mask','released','technical_QC_pass','policy_prediction_available'):
        np.testing.assert_array_equal(cpu[name],gpu[name])

def test_prediction_preserves_original_cuda_error(fitted,monkeypatch):
    import copy
    learner=copy.deepcopy(fitted['b']['models'][fitted['b']['primary']]);learner.device='cuda'
    def fail(*args,**kwargs):raise RuntimeError('intentional CUDA kernel error')
    monkeypatch.setattr(learner,'predict_batch',fail)
    with pytest.raises(RuntimeError,match='intentional CUDA kernel error'):
        learner.predict(fitted['x'][:2],fitted['p'].iloc[:2])
    assert learner.device=='cuda' and next(learner.network.parameters()).device.type=='cpu'

def test_real_cuda_3000_assays_and_peak(tmp_path):
    torch=m.torch
    if not torch.cuda.is_available():pytest.fail('CUDA is required on this acceptance host')
    torch.cuda.reset_peak_memory_stats()
    start=time.monotonic();m.s7_model_preflight(config('--device','cuda'),features=3000)
    record=dict(features=3000,actual_device='cuda',forward_backward=True,wall_seconds=time.monotonic()-start,
        peak_allocated_bytes=torch.cuda.max_memory_allocated(),peak_reserved_bytes=torch.cuda.max_memory_reserved(),gpu=torch.cuda.get_device_name())
    (tmp_path/'gpu-3000-evidence.json').write_text(json.dumps(record,indent=2))


if __name__ == '__main__':
    if sys.argv[1:2] == ['--smoke-data']:
        make_smoke_data(sys.argv[2:])
    else:
        raise SystemExit(pytest.main([__file__, *sys.argv[1:]]))
