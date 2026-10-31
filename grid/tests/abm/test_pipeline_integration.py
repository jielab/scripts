"""Actual adapters on isolated synthetic input files, without LD or benchmark fits."""
import json
import shutil
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import joblib
import numpy as np
import pandas as pd
from abm_core import pipeline
from abm_core.derived_report import report_from
from abm_core.deploy import assess_dataset
from abm_core.model import file_sha256,write_json
from abm_core.contracts import write_array_manifest
from abm_core.representation import BlockGeometry
from browser_index import Index


def prepared_fixture(root):
    rng=np.random.default_rng(20261020);n=3600
    features=['csx.'+p for p in ('AFR','EAS','EUR','SAS')]+[f'match_PC{i}' for i in range(1,11)]+['age','sex','PC1','PC2']
    x=rng.normal(size=(n,len(features)))
    d=pd.DataFrame(x,columns=features)
    d['eid']=[f'{i:06d}' for i in range(n)];d['family_id']=['g'+str(i//2) for i in range(n)]
    d['split']=np.where(np.arange(n)<3200,'train','test');d['ancestry']='EUR'
    y=2*x[:,0]+.4*x[:,14]+rng.normal(size=n)
    groups=dict(csx=features[:4],ancestry=features[4:14],covariates=features[14:])
    tables={t:d.assign(y=((y>0).astype(float) if t=='t2dm' else y*(2 if t=='ldl' else 1))) for t in ('height','ldl','t2dm')}
    for t in tables:tables[t].loc[n-5:,'y']=np.nan
    path=root/'input/_LOCAL_ONLY/prepared.joblib';path.parent.mkdir(parents=True)
    joblib.dump(dict(state='features_ready',tables=tables,cohort=d[['eid','family_id','split']],
        feature_groups={t:groups for t in tables},metadata=dict(grouping_status='synthetic_family_groups',synthetic_only=True)),path)
    config=json.loads((pipeline.ROOT/'config/abm.json').read_text())
    config['block_config']=str(root/'NONEXISTENT_LD_CONFIG.json')
    config['head_options'].update(n_bins=2,n_boot=30,min_calibration_groups=2)
    config['evaluation'].update(bootstrap_replicates=30,min_matched_pairs=2,min_matched_groups_each_side=2)
    config['support'].update(radius_k_ref=20,min_count=2,min_group_ess=2,query_batch=128)
    config['predictors'].update(moe_clusters=[4],ridge_grid=[10.],logistic_c_grid=[.1])
    p=root/'abm.json';p.write_text(json.dumps(config))
    return path,p


class RealAdapterTests(unittest.TestCase):
    def test_binary_assessment_loss_change_reuses_scientific_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);prepared,cfg=prepared_fixture(root)
            args=['--analysis-root',td,'--prepared',str(prepared),'--config',str(cfg),'--traits','t2dm',
                  '--predictor','moe','--representation','csx_scores','--device','cpu']
            pipeline.main(args)
            p=root/'abm/FINAL/t2dm'
            first=p/json.loads((p/'moe.CURRENT.json').read_text())['run']
            config=json.loads(cfg.read_text());config['primary_loss']['t2dm']='logloss';cfg.write_text(json.dumps(config))
            with patch.object(pipeline,'fit_predictors',side_effect=AssertionError('loss does not retrain fixed predictors')), \
                 patch.object(pipeline,'count_twins',side_effect=AssertionError('loss does not rebuild label-free support')):
                pipeline.main(args)
            second=p/json.loads((p/'moe.CURRENT.json').read_text())['run']
            from abm_core.model import ABMAssessor
            self.assertEqual(ABMAssessor.load(second/'_LOCAL_ONLY/assessor.json').heads.loss_name,'logloss')
            self.assertNotEqual(first,second)
            self.assertEqual(file_sha256(first/'_LOCAL_ONLY/predictor_bundle.joblib'),
                             file_sha256(second/'_LOCAL_ONLY/predictor_bundle.joblib'))
            cal=pd.read_csv(second/'test_count_calibration.tsv',sep='\t')
            self.assertTrue(cal.assessment_loss.eq('logloss').all())

    def test_blocknet_actual_producer_to_seven_stages_cpu(self):
        import sys
        sys.path.insert(0,str(pipeline.ROOT/'tests'))
        from test_block_producer import fixture as block_fixture
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);prepared,cfg=prepared_fixture(root)
            bdir=root/'block_input';bdir.mkdir()
            block_cfg,_=block_fixture(bdir,'ALL.37.bed',n=3600)
            import block_inputs
            block_inputs.prepare_trait(block_cfg,'height',device='cpu')
            config=json.loads(cfg.read_text());config['block_config']=str(block_cfg);cfg.write_text(json.dumps(config))
            pipeline.main(['--analysis-root',td,'--prepared',str(prepared),'--config',str(cfg),'--traits','height',
                           '--device','cpu','--smoke','--smoke-rows','160'])
            parent=root/'abm/FINAL-SMOKE/height'
            self.assertFalse((parent/'blocknet.CURRENT.json').exists())
            run=next(p.parent for p in parent.glob('*/COMPLETE.json'))
            self.assertEqual(set(json.loads((run/'PHASES.json').read_text())['phases']),set(pipeline.PHASES))
            self.assertTrue(json.loads((run/'COMPLETE.json').read_text())['smoke'])
            self.assertIn(json.loads((run/'block_staging.json').read_text())['strategy'],['ram','mmap'])
            from browser_index import Index
            Index(root,root/'final').run(['height'])
            catalogue=pd.read_csv(root/'final/catalogue.tsv',sep='\t')
            self.assertFalse(catalogue.module.eq('abm').any())
            self.assertFalse(catalogue.source.str.contains('SMOKE').any())

    def test_seven_stages_score_only_loss_cache_deploy_report_and_shiny(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);prepared,cfg=prepared_fixture(root)
            args=['--analysis-root',td,'--prepared',str(prepared),'--config',str(cfg),
                  '--traits','height','--predictor','moe','--representation','csx_scores','--device','cpu']
            # Full production adapters execute; only unused dependencies are guarded.
            with patch.object(pipeline,'locate_blocks',side_effect=AssertionError('LD must not be read')), \
                 patch('prsformer_resume.dispatch',side_effect=AssertionError('benchmark must not train')):
                pipeline.main(args)
            parent=root/'abm/FINAL/height';run=parent/json.loads((parent/'moe.CURRENT.json').read_text())['run']
            state=json.loads((run/'PHASES.json').read_text())
            self.assertEqual(set(state['phases']),set(pipeline.PHASES))
            self.assertEqual(len(pd.read_parquet(run/'_LOCAL_ONLY/test/predictions.parquet')),400)
            self.assertEqual(len(pd.read_csv(run/'count_incremental_contrasts.tsv',sep='\t')),4)
            original={str(p.relative_to(run)):file_sha256(p) for p in run.rglob('*') if p.is_file()}
            with patch.object(pipeline,'fit_predictors',side_effect=AssertionError('report must not fit')):
                report=report_from(run,root/'abm_reports/derived')
            self.assertEqual(original,{str(p.relative_to(run)):file_sha256(p) for p in run.rglob('*') if p.is_file()})
            self.assertEqual(len(pd.read_csv(report/'count_incremental_contrasts.tsv',sep='\t')),4)
            # A run that never exported diagnostics may fit those development-only
            # diagnostics; neither primary selected masks nor predictors change.
            incomplete=root/'without_diagnostic_export';shutil.copytree(run,incomplete)
            (incomplete/'_LOCAL_ONLY/diagnostics.joblib').unlink()
            phases=json.loads((incomplete/'PHASES.json').read_text())
            phases['phases']['fit_assessor'].pop('_LOCAL_ONLY/diagnostics.joblib')
            write_json(incomplete/'PHASES.json',phases)
            with patch.object(pipeline,'fit_predictors',side_effect=AssertionError('must not fit predictors')):
                derived=report_from(incomplete,root/'abm_reports/development_diagnostics')
            added=pd.read_parquet(derived/'_LOCAL_ONLY/test/predictions.parquet')
            frozen=pd.read_parquet(run/'_LOCAL_ONLY/test/predictions.parquet')
            for column in ['selected','switched','p_policy']:
                pd.testing.assert_series_equal(added[column],frozen[column])
            # ID-aligned score-only deployment reproduces frozen support, also after query reordering.
            features=pd.read_parquet(run/'_LOCAL_ONLY/features.parquet')
            query=features[features.role.eq('test')].reset_index(drop=True)
            bundle=joblib.load(run/'_LOCAL_ONLY/predictor_bundle.joblib')
            geometry=json.loads((run/'geometry_manifest.json').read_text())
            transform=BlockGeometry(**geometry['state'])
            results=[]
            for i,frame in enumerate((query,query.sample(frac=1,random_state=3).reset_index(drop=True))):
                z=transform.transform(frame[[bundle['features'][j] for j in bundle['csx_columns']]].to_numpy(float))
                arr=root/f'geometry{i}.npy';np.save(arr,z)
                binding=write_array_manifest(arr,frame,geometry_sha256=geometry['geometry_sha256'])
                out,_=assess_dataset(run,frame,arr,binding,device='cpu');results.append(out.set_index('id').sort_index())
            for col in ('n_twins','supported','selected','p_policy'):
                np.testing.assert_allclose(results[0][col],results[1][col],rtol=1e-12,atol=1e-12)
            # Changing only test Y creates a truthful new run but reuses predictor and support caches.
            saved=joblib.load(prepared);saved['tables']['height'].loc[lambda d:d.split.eq('test'),'y']=5.
            joblib.dump(saved,prepared)
            with patch.object(pipeline,'fit_predictors',side_effect=AssertionError('test Y must not refit predictor')), \
                 patch.object(pipeline,'count_twins',side_effect=AssertionError('test Y must not recompute support')):
                pipeline.main(args)
            run2=parent/json.loads((parent/'moe.CURRENT.json').read_text())['run']
            a=pd.read_parquet(run/'_LOCAL_ONLY/support.parquet');b=pd.read_parquet(run2/'_LOCAL_ONLY/support.parquet')
            pd.testing.assert_frame_equal(a,b)
            # A real score-value change invalidates both fit and geometric support.
            saved['tables']['height'].loc[0,'csx.AFR']+=.3;joblib.dump(saved,prepared)
            with patch.object(pipeline,'fit_predictors',wraps=pipeline.fit_predictors) as fit, \
                 patch.object(pipeline,'count_twins',wraps=pipeline.count_twins) as scan:
                pipeline.main(args)
            self.assertEqual(fit.call_count,1);self.assertEqual(scan.call_count,1)
            from final_compare import run as final_summary
            final_summary(root,['height'])
            combined=pd.read_csv(root/'final/count_incremental_contrasts.tsv',sep='\t')
            self.assertEqual(set(combined.report_scope),{'derived_report'})
            self.assertEqual(len(pd.read_csv(root/'final/performance.tsv',sep='\t')),0)
            self.assertGreater(len(pd.read_csv(root/'final/sensitivity_performance.tsv',sep='\t')),0)
            # Actual R/Shiny runtime exercises the production index, plot and reactive page.
            db=Index(root,root/'final').run(['height'])['db']
            from twin import r_runtime
            try:r=r_runtime()
            except RuntimeError:self.skipTest('R/Shiny runtime absent: actual rendering not executed')
            result=subprocess.run([*r,'--vanilla',str(pipeline.ROOT/'tests/shiny_fixture.R'),db,
                str(pipeline.ROOT/'shiny')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+'\n'+result.stderr)
            # Corruption is rejected before deriving any report or reading truth.
            with (run/'_LOCAL_ONLY/truth.parquet').open('ab') as f:f.write(b'tamper')
            with self.assertRaisesRegex(ValueError,'artifact changed'):report_from(run,root/'abm_reports/reject')

if __name__=='__main__':unittest.main()
