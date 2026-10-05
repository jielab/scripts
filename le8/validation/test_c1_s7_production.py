"""S7 production contracts T01–T52; real fits use only synthetic /tmp data."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
if 'c1_abm' in sys.modules:
	m=sys.modules['c1_abm']
else:
	spec=importlib.util.spec_from_file_location('c1_abm',ROOT/'f/c1.abm.py')
	m=importlib.util.module_from_spec(spec);sys.modules['c1_abm']=m;spec.loader.exec_module(m)


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
	report=source.replace("title='Frozen S7 models: identical held-out participants'", "title='Updated figure caption'")
	assert source!=report
	assert m.s7_numerical_code_hash(source)==m.s7_numerical_code_hash(report)
	assert m.s7_numerical_code_hash(source)!=m.s7_numerical_code_hash(source.replace('a.s7_reconstruction_weight = .05','a.s7_reconstruction_weight = .07'))
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
	old=m.configure(m.reference_parser().parse_args([]))
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
