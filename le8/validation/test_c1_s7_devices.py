"""T21/T24: actual reconstruction-free CPU/CUDA optimization and CPU restore."""
import numpy as np
import pandas as pd
import pytest
from threadpoolctl import threadpool_limits

from test_c1_s7_production import m, config


@pytest.mark.parametrize('device',['cpu','cuda'])
def test_real_supervised_qk_gradients_and_cpu_restore(device):
	if device=='cuda' and not m.torch.cuda.is_available():
		pytest.skip('Real CUDA unavailable; never substitute CPU and call it CUDA')
	a=config('--device',device,'--epochs','2','--batch-size','128')
	a.s7_reconstruction_weight=0.;a.group_col=''
	rng=np.random.default_rng(4);n=500
	x=rng.normal(size=(n,8));p=pd.DataFrame(dict(eid=[f'device-{i}' for i in range(n)],age=rng.uniform(40,70,n),sex=rng.integers(0,2,n)))
	time=rng.exponential(8,n);event=(time<9).astype(int);time=np.minimum(time,9)
	data=m.S7FitData(x,p,time,event,m.s7_schema([f'F{i}' for i in range(8)],a),'build')
	with threadpool_limits(2):
		m.torch.set_num_threads(2)
		learner=m.S7AttentionLearner().fit(data.take(np.arange(400)),data.take(np.arange(400,500),'tune_model'),np.ones(400),vars(a))
		assert {'query.weight','key.weight'}<=learner.gradient_names
		assert any(name.startswith('context.') for name in learner.gradient_names)
		assert learner.fit_status['selected_tensor_changed_count']>0
		assert learner.fit_status['supervised_steps']>=2
		assert learner.fit_status['actual_device']==device
		restored=m.S7AttentionLearner.restore_state(learner.export_state(),'cpu')
		full=learner.predict(x[400:],p.iloc[400:]).probability
		one=restored.predict(x[400:401],p.iloc[400:401]).probability
		np.testing.assert_allclose(one,full[:1],atol=1e-6,rtol=0)
	print(dict(device=device,**learner.fit_status))
