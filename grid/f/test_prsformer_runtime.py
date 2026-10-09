"""Opt-in tests of native CUDA attention and the complete three-trait wrapper.

GRID_TEST_PRSFORMER_GPU=1 ~/.venvs/grid-prsformer/bin/python -m unittest discover \
    -s f -p 'test_prsformer_runtime.py' -v
Set GRID_PRSFORMER_TEST_WORK to retain synthetic fixtures/logs for inspection.
These tests never read UKB inputs or use the production cache/output directories.
"""
from pathlib import Path
import importlib.util
import json
import os
import signal
import subprocess
import time
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prsformer_gpu_test', ROOT/'f/3.prsformer.py')
prs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prs)
ENABLED = os.environ.get('GRID_TEST_PRSFORMER_GPU') == '1'
UPSTREAM = Path(os.environ.get('PRSFORMER_UPSTREAM_DIR', '/mnt/f/software/PRSformer'))


@unittest.skipUnless(ENABLED, 'Set GRID_TEST_PRSFORMER_GPU=1 to exercise the real CUDA kernels')
class NativeAttention(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import torch
        cls.torch = torch
        if not torch.cuda.is_available():
            raise RuntimeError('GPU tests requested but CUDA unavailable')
        torch.set_num_threads(4)
        cls.natten, version = prs.neighborhood_runtime()
        if version != '0.21.7':
            raise RuntimeError('These compatibility tests require the documented modern NATTEN stack')

    def test_fused_outputs_and_gradients_match_independent_sliding_attention(self):
        torch = self.torch
        for length, kernel, dilation in [(13,3,1),(17,3,2),(19,5,2),(389,385,1)]:
            for dtype in (torch.float32, torch.float16, torch.bfloat16):
                with self.subTest(length=length, kernel=kernel, dilation=dilation, dtype=dtype):
                    torch.manual_seed(812)
                    layer = prs.modern_natten_compatibility(self.natten).NeighborhoodAttention1D(
                        dim=64, num_heads=4, kernel_size=kernel, dilation=dilation).cuda().to(dtype)
                    x = torch.randn(2,length,64,device='cuda',dtype=dtype,requires_grad=True)
                    actual = layer(x)
                    # Independently construct the shifted full window within each
                    # dilation residue class, including left/right edge windows.
                    allowed = torch.zeros(length,length,dtype=torch.bool,device='cuda')
                    for i in range(length):
                        neighbors = list(range(i % dilation,length,dilation))
                        start = max(0,min(neighbors.index(i)-kernel//2,len(neighbors)-kernel))
                        allowed[i,neighbors[start:start+kernel]] = True
                    q,k,v = layer.qkv(x).reshape(2,length,3,4,16).unbind(2)
                    logits = torch.einsum('bihd,bjhd->bhij',q.float(),k.float()) * .25
                    weights = logits.masked_fill(~allowed,-torch.inf).softmax(-1)
                    attended = torch.einsum('bhij,bjhd->bihd',weights,v.float()).reshape(2,length,64)
                    expected = layer.proj(attended.to(dtype))
                    tolerance = 3e-5 if dtype == torch.float32 else (.004 if dtype == torch.float16 else .035)
                    torch.testing.assert_close(actual,expected,atol=tolerance,rtol=tolerance)
                    probe = torch.randn_like(actual)
                    variables = (x,*layer.parameters())
                    grad_actual = torch.autograd.grad((actual*probe).float().mean(),variables,retain_graph=True)
                    grad_expected = torch.autograd.grad((expected*probe).float().mean(),variables)
                    for a,b in zip(grad_actual,grad_expected):
                        self.assertTrue(bool(torch.isfinite(a).all()))
                        torch.testing.assert_close(a,b,atol=tolerance/10,rtol=tolerance*2)
                    self.assertEqual(set(layer.state_dict()), {'qkv.weight','qkv.bias','proj.weight','proj.bias'})

    def test_unsupported_old_options_fail(self):
        adapter = prs.modern_natten_compatibility(self.natten)
        for options in ({'rel_pos_bias':True},{'attn_drop':.1}):
            with self.assertRaisesRegex(ValueError,'rel_pos_bias=False'):
                adapter.NeighborhoodAttention1D(dim=64,num_heads=4,kernel_size=3,**options)
        with self.assertRaisesRegex(ValueError,'always requires'):
            adapter.use_fused_na(False)

    def test_default_runtime_without_input_files_or_output_directory(self):
        with tempfile.TemporaryDirectory(prefix='prsformer-runtime-test-') as td:
            work=Path(td)
            command=['bash',str(ROOT/'3.prsformer.sh'),'--python',sys.executable,
                     '--upstream-dir',str(UPSTREAM),'--check-runtime',
                     '--cache-dir',str(work/'absent-cache'),'--run-dir',str(work/'absent-run')]
            result=subprocess.run(command,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('official model forward/backward passed',result.stdout)
            self.assertFalse((work/'absent-cache').exists())
            self.assertFalse((work/'absent-run').exists())


@unittest.skipUnless(ENABLED, 'Set GRID_TEST_PRSFORMER_GPU=1 to run native GPU integration')
class ThreeTraitPipeline(unittest.TestCase):
    def test_prepare_train_predict_resume_and_publish(self):
        import pyreadr
        persistent=os.environ.get('GRID_PRSFORMER_TEST_WORK')
        if persistent:
            work=Path(persistent).resolve(); work.mkdir(parents=True,exist_ok=True)
        else:
            temporary=tempfile.TemporaryDirectory(prefix='prsformer-pipeline-test-')
            self.addCleanup(temporary.cleanup); work=Path(temporary.name)
        plink=os.environ.get('GRID_TEST_PLINK2',str(Path.home()/'miniforge3/envs/grid/bin/plink2'))
        n,m=160,1024
        rng=np.random.default_rng(20261008)
        ids=[str(200000+i) for i in range(n)]
        g=rng.binomial(2,.32,(n,m)).astype(np.int8)
        g[3,5]=-1
        g[:,0]=0  # A monomorphic SNP must be removed by train-only QC.
        (work/'gen').mkdir(exist_ok=True)
        for chrom in (1,2):
            prefix=work/'gen'/f'bed{chrom}'
            order=np.arange(n) if chrom==1 else rng.permutation(n)
            variants=range((chrom-1)*512,chrom*512)
            prefix.with_suffix('.fam').write_text(''.join(f'0 {ids[i]} 0 0 0 -9\n' for i in order))
            prefix.with_suffix('.bim').write_text(''.join(f'{chrom} rs{j} 0 {j+100} C A\n' for j in variants))
            packed=bytearray(b'\x6c\x1b\x01'); bits={-1:1,0:3,1:2,2:0}
            for j in variants:
                row=g[order,j]
                for start in range(0,n,4):
                    packed.append(sum(bits[int(v)] << (2*k) for k,v in enumerate(row[start:start+4])))
            prefix.with_suffix('.bed').write_bytes(packed)
            subprocess.run([plink,'--bfile',str(prefix),'--make-pgen','--out',str(work/'gen'/f'chr{chrom}')],
                           check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        age=rng.normal(50,5,n)
        pheno=pd.DataFrame({'eid':ids,'age':age,'sex':rng.integers(0,2,n),
            'PC1':rng.normal(size=n),'PC2':rng.normal(size=n),
            'genetic_ancestry':np.where(np.arange(n)%2,'AFR','EUR'),
            'height':165+.2*age+.2*g[:,1]+rng.normal(size=n),
            'ldl':3+.1*g[:,2]+rng.normal(0,.3,n),
            't2dm.Yr2e':(np.arange(n)%3==0).astype(float),'t2dm.Yt2e':0.})
        pheno.loc[[4,10,55,145],'ldl']=np.nan
        pyreadr.write_rds(str(work/'pheno.rds'),pheno,compress='gzip')
        pd.DataFrame({'SNP':[f'rs{j}' for j in range(m)],'CHR':np.repeat([1,2],512),
                      'BP':np.arange(m)+100}).to_csv(work/'snps.tsv',sep='\t',index=False)
        base=['bash',str(ROOT/'3.prsformer.sh'),'--python',sys.executable,
              '--upstream-dir',str(UPSTREAM),'--traits','height,ldl,t2dm',
              '--dir-gen',str(work/'gen'),'--chrs','1-2','--pheno-file',str(work/'pheno.rds'),
              '--ancestry-file','none','--snp-list',str(work/'snps.tsv'),'--remove','none',
              '--cache-dir',str(work/'cache'),'--output-root',str(work/'reports'),
              '--score-dir',str(work/'scores'),'--epochs','2','--patience','2']
        def run(label, extra=(), success=True):
            result=subprocess.run(base+list(extra),text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            (work/f'{label}.log').write_text(result.stdout)
            if success: self.assertEqual(result.returncode,0,f'{label}: {result.stdout[-10000:]}')
            else: self.assertNotEqual(result.returncode,0)
            return result
        run('all')
        meta=json.loads((work/'cache/prepare.json').read_text())
        self.assertEqual(meta['shape'],[n,m-1])
        np.testing.assert_array_equal(np.load(work/'cache/genotypes.npy'),g[:,1:])
        original=pd.read_csv(work/'cache/run/test_predictions.tsv.gz',sep='\t')
        run('check',['--check'])
        refused=run('overwrite-rejected',success=False)
        self.assertIn('saved model already exists',refused.stdout)
        completed=run('completed-resume-rejected',['--mode','train','--resume'],success=False)
        self.assertIn('needs both training.pt and model.pt',completed.stdout)
        run('predict',['--mode','predict','--checkpoint',str(work/'cache/run/model.pt'),
                       '--run-dir',str(work/'predict')])
        predicted=pd.read_csv(work/'predict/test_predictions.tsv.gz',sep='\t')
        pd.testing.assert_frame_equal(original,predicted,check_exact=False,rtol=1e-6,atol=1e-7)
        run('report',['--mode','report','--replace'])
        # Interrupt only the isolated child process group after an atomic epoch
        # checkpoint exists, then use the normal shell entry point to resume.
        resume_args=['--mode','train','--run-dir',str(work/'resume-run'),'--epochs','4','--patience','4']
        with (work/'interrupted.log').open('w') as log:
            proc=subprocess.Popen(base+resume_args,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            deadline=time.monotonic()+90
            try:
                while not (work/'resume-run/training.pt').exists():
                    if proc.poll() is not None or time.monotonic()>deadline:
                        self.fail('Interrupted-run fixture did not reach an epoch checkpoint: '+(work/'interrupted.log').read_text())
                    time.sleep(.02)
            finally:
                if proc.poll() is None:
                    os.killpg(proc.pid,signal.SIGTERM)
                proc.wait(timeout=15)
        interrupted_state=prs.import_torch().load(work/'resume-run/training.pt',map_location='cpu',weights_only=True)
        run('resume',resume_args+['--resume'])
        restored_history=pd.read_csv(work/'resume-run/history.tsv',sep='\t')
        old_history=pd.DataFrame(interrupted_state['history'])
        np.testing.assert_allclose(restored_history.loss.iloc[:len(old_history)],old_history.loss,rtol=1e-12)
        self.assertEqual(restored_history.epoch.max(),4)
        run('uninterrupted',['--mode','train','--run-dir',str(work/'uninterrupted'),
                             '--epochs','4','--patience','4'])
        resumed=pd.read_csv(work/'resume-run/test_predictions.tsv.gz',sep='\t')
        continuous=pd.read_csv(work/'uninterrupted/test_predictions.tsv.gz',sep='\t')
        # Fused CUDA backward uses parallel reductions; independent FP16 runs
        # need not be bitwise identical, even before the interruption.
        exact_columns=[c for c in resumed if c not in ('prediction','genetic_score')]
        pd.testing.assert_frame_equal(resumed[exact_columns],continuous[exact_columns])
        for column in ('prediction','genetic_score'):
            np.testing.assert_allclose(resumed[column],continuous[column],rtol=2e-3,atol=3e-3)
        for trait in ('height','ldl','t2dm'):
            score=next(iter(pyreadr.read_r(str(work/'scores'/trait/'3.prsformer.scores.rds')).values()))
            self.assertGreater(len(score),0)
        for kind in ('performance','training'):
            for extension in ('xlsx','png'):
                self.assertTrue((work/'reports/prsformer'/f'3.prsformer.{kind}.{extension}').is_file())
        self.assertEqual(subprocess.run(['git','-C',str(UPSTREAM),'diff','--exit-code'],capture_output=True).returncode,0)
        print(f'GPU three-trait pipeline artifacts: {work}',flush=True)


if __name__ == '__main__':
    unittest.main()
