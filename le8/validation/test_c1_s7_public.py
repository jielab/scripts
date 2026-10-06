"""T30/T45/T48–T50: actual le8.sh transactions, cache, projection and Final."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import shlex
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def test_public_cpu_frozen_roundtrip_and_final():
	root=Path(os.getenv('LE8_S7_PUBLIC_DIR') or tempfile.mkdtemp(prefix='s7-public-',dir='/tmp'))
	inputs=Path(os.getenv('LE8_S7_INPUT_DIR') or root.parent/(root.name+'-synthetic-input'))
	root.mkdir(parents=True,exist_ok=True)
	if not (inputs/'phe.csv').exists():
		subprocess.run([sys.executable,str(ROOT/'validation/make_s7_smoke_data.py'),'--out',str(inputs)],check=True)
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
	assert (out/'test_individuals.rds').is_file()
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
	from test_c1_s7_production import m, config, synthetic
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
		private_RDS_roundtrip=True,outcome_free_external=True,training_inputs_unavailable=True,model_sha256=sha(),
		final_registered_models=int(registry.model_id.nunique()),synthetic_only=True),indent=2))
