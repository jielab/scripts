"""Scheduling tests use fake bsub; no analysis jobs are submitted."""
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('hpc_batch', Path(__file__).resolve().parents[1] / 'hpc_batch.py')
batch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(batch)


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.psam = self.root / 'samples.psam'
        self.psam.write_text('#FID IID SEX\n' + ''.join(f'0 s{i} {1+i%2}\n' for i in range(25)))
        (self.root / 'f').mkdir()
        (self.root / 'f/hpc.runtime.env').write_text('# test runtime\n')
        self.args = ['ibdmix', '--target', 'ukb', '--dir-out', str(self.root / 'ukb'),
                     '--postfix', '[001-003]', '--keep-psam', f'{self.psam},10,[001-003]', '--keep-males']

    def tearDown(self):
        self.tmp.cleanup()

    def test_default_result_locations_follow_analysis_root(self):
        with patch.dict(os.environ, {'GU_REF_ROOT': '/work/sph-huangj/data/gen', 'GU_ANALYSIS_BASE': '/work/sph-huangj/analysis/gu'}):
            self.assertEqual(str(batch.result_base(['ibdmix','--target','ukb'])), '/work/sph-huangj/analysis/gu/ukb/37')
            self.assertEqual(str(batch.result_base(['ibdmix','--dir-target','/mnt/f/gen/1kg/37/pfile/chr'])), '/work/sph-huangj/analysis/gu/1kg/37')
            self.assertEqual(str(batch.result_base(['ibdmix','--target','ukb','--postfix','013'])), '/work/sph-huangj/analysis/gu/ukb/37/cohort')

    def test_missing_rscript_fails_before_submission(self):
        with patch.dict(os.environ, {'GU_RESULTS_RSCRIPT':str(self.root/'missing/Rscript')}):
            with self.assertRaisesRegex(ValueError, 'missing .*/missing/Rscript'):
                batch.preflight_publication()

    def test_running_snapshot_keeps_helper_after_live_edit(self):
        project=self.root/'scripts/gu';(project/'f').mkdir(parents=True)
        helper=project/'f/ukb_prepare.py';helper.write_text('original helper')
        plan=self.root/'plan';plan.mkdir()
        with patch.object(batch,'ROOT',project): frozen=batch.snapshot_code(plan)
        helper.unlink()
        self.assertEqual((frozen/'f/ukb_prepare.py').read_text(),'original helper')
        self.assertTrue((plan/'code.sha256.json').is_file())

    def test_padding_and_lists(self):
        self.assertEqual(batch.expand('[013-051]'), [f'{i:03}' for i in range(13, 52)])
        self.assertEqual(batch.expand('013,015,049'), ['013', '015', '049'])

    def test_invalid_selectors(self):
        for text in ('[051-013]', '[0-2]', '[013-015', '013,013', '../013'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                batch.expand(text)

    def test_chunk_identity_and_male_flag(self):
        pairs = batch.cohorts(self.args)
        self.assertEqual([p for p, _ in pairs], ['001', '002', '003'])
        for n, (_, cmd) in enumerate(pairs, 1):
            self.assertEqual(batch.get(cmd, '--keep-psam'), f'{self.psam},10,{n}')
            self.assertIn('--keep-males', cmd)

    def test_mismatch_and_empty_chunks(self):
        for flag, value in [('--postfix', '[001-002]'), ('--postfix', '[002-004]'),
                            ('--keep-psam', f'{self.psam},10,[003-005]')]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                batch.cohorts(batch.put(self.args, flag, value))

    def test_single_non_numeric_postfix(self):
        args = batch.put(batch.put(self.args, '--postfix', 'compare'), '--keep-psam', f'{self.psam},10,2')
        self.assertEqual(batch.cohorts(args)[0][0], 'compare')

    def test_disjoint_list(self):
        args = batch.put(batch.put(self.args, '--postfix', '001,003'), '--keep-psam', f'{self.psam},10,[001,003]')
        self.assertEqual([p for p, _ in batch.cohorts(args)], ['001', '003'])

    def run_main(self, args, fake):
        with patch.object(batch, 'ROOT', self.root), patch.dict(os.environ, {'GU_ANALYSIS_BASE': str(self.root)}), patch.object(batch.subprocess, 'run', fake), patch.object(batch, 'preflight_axt'), patch.object(batch, 'preflight_inputs'), patch.object(batch, 'preflight_publication'), patch.object(batch, 'snapshot_code', return_value=self.root):
            return batch.main(args)

    def test_dry_run_never_submits_and_worker_runs(self):
        with patch.object(batch.subprocess, 'run') as fake:
            self.run_main(self.args + ['--run-cmd', 'FALSE'], fake)
            fake.assert_not_called()
        files = sorted(self.root.glob('ukb.jobs/*/*.cmd'))
        self.assertEqual(len(files), 3)
        text = files[0].read_text()
        self.assertIn('--run-cmd TRUE', text)
        self.assertIn('--foreground TRUE', text)
        self.assertIn('--GPU FALSE', text)
        self.assertIn('${LSB_JOBID:', text)
        self.assertIn('--postfix 001', text)

    def test_submits_independently(self):
        from unittest.mock import Mock
        fake = Mock(side_effect=[subprocess.CompletedProcess([], 0, f'Job <{n}> is submitted', '') for n in (11, 12, 13)])
        self.run_main(self.args, fake)
        self.assertEqual(fake.call_count, 3)
        self.assertTrue(all(call.args == (['bsub'],) for call in fake.call_args_list))

    def test_explicit_gpu_plan_requests_one_gpu_without_probing_login_node(self):
        with patch.dict(os.environ, {'submit':'hpc', 'GU_GPU_QUEUE':'4a100-40', 'GU_GPU_MEMORY_GIB':'160'}), patch.object(batch.subprocess, 'run') as fake:
            self.run_main(self.args + ['--GPU', 'TRUE', '--run-cmd', 'FALSE'], fake)
            fake.assert_not_called()
        text = next(self.root.glob('ukb.jobs/*/*.cmd')).read_text()
        self.assertIn('#BSUB -gpu "num=1/host:mode=exclusive_process:j_exclusive=yes"', text)
        self.assertIn('#BSUB -q 4a100-40', text)
        self.assertIn('#BSUB -n 16', text)
        self.assertIn('#BSUB -M 160GB', text)
        self.assertIn('--GPU TRUE', text)

    def test_explicit_cpu_overrides_site_gpu_default(self):
        with patch.dict(os.environ, {'submit':'hpc', 'GU_IBDMIX_GPU_DEFAULT':'TRUE'}), patch.object(batch.subprocess, 'run') as fake:
            self.run_main(self.args + ['--GPU', 'FALSE', '--run-cmd', 'FALSE'], fake)
        text = next(self.root.glob('ukb.jobs/*/*.cmd')).read_text()
        self.assertNotIn('#BSUB -gpu', text)
        self.assertIn('--GPU FALSE', text)

    def test_submission_failure_stops_and_records(self):
        from unittest.mock import Mock
        fake = Mock(side_effect=[subprocess.CompletedProcess([], 0, 'Job <11> is submitted', ''), subprocess.CompletedProcess([], 1, '', 'queue unavailable')])
        with self.assertRaises(RuntimeError):
            self.run_main(self.args, fake)
        self.assertEqual(fake.call_count, 2)
        journal = next(self.root.glob('ukb.jobs/*/submission-*.jsonl')).read_text()
        self.assertIn('"job_id": "11"', journal)
        self.assertIn('queue unavailable', journal)


if __name__ == '__main__':
    unittest.main()
