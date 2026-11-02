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
        with patch.object(batch, 'ROOT', self.root), patch.dict(os.environ, {'GU_ANALYSIS_BASE': str(self.root)}), patch.object(batch.subprocess, 'run', fake), patch.object(batch, 'preflight_axt'):
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
