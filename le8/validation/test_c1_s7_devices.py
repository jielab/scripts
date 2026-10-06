"""T21/T24: actual reconstruction-free CPU/CUDA optimization and CPU restore."""
import numpy as np
import pandas as pd
import pytest
from threadpoolctl import threadpool_limits
from unittest import mock

from test_c1_s7_production import m, config


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
