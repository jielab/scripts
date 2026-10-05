"""Five actual learner experiments; set LE8_S7_SIMULATION_DIR to retain artifacts."""
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
import pytest
from threadpoolctl import threadpool_limits

from test_c1_s7_production import m, config


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
