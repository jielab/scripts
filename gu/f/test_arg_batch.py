"""Scheduling tests: all scheduler calls are mocked; no scientific jobs run."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import arg_batch as batch


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / 'gen/1kg/37').mkdir(parents=True)
        self.env = patch.dict(os.environ, {'GU_REF_ROOT': str(self.root / 'gen'),
                                          'GU_ANALYSIS_BASE': str(self.root / 'analysis'),
                                          'ARG_MEMORY_LIMIT_GB': '32', 'GU_QUEUE': '33'})
        self.env.start()
        self.output = contextlib.redirect_stdout(io.StringIO())
        self.output.__enter__()

    def tearDown(self):
        self.output.__exit__(None, None, None)
        self.env.stop()
        self.temp.cleanup()

    def snapshot(self, run):
        frozen = run / 'code/gu'
        (frozen / 'f').mkdir(parents=True)
        return frozen

    def make_plan(self, extra=(), submit=False):
        args = ['build', '--grch', '37', '--target', '1kg', '--chr', '1,20,X',
                '--method', 'threads,needle,tsinfer,singer', '--format', 'trace',
                '--threads', '8', '--jobs', '1', '--run-dir', str(self.root / 'run')]
        if not submit:
            args.append('--plan-only')
        with patch.object(batch, 'snapshot_code', side_effect=self.snapshot), patch.object(batch.submit_jobs, 'run') as runner:
            run = batch.plan(args + list(extra))
            if submit:
                self.assertEqual(len(runner.call_args.args[0]), 12)
            else:
                runner.assert_not_called()
        return run, json.loads((run / 'plan.json').read_text())

    def test_matrix_and_worker_scripts(self):
        run, plan = self.make_plan()
        tasks = plan['tasks']
        self.assertEqual({(t['method'], t['chromosome']) for t in tasks},
                         {(m, c) for m in batch.METHODS for c in ('1', '20', 'X')})
        for field in ('output', 'prepared', 'directory', 'cmd'):
            self.assertEqual(len({t[field] for t in tasks}), 12)
        self.assertEqual(plan['slots'], 10)  # existing QiMing quarter-node rule
        for t in tasks:
            self.assertEqual(batch.get(t['argv'], '--method'), t['method'])
            self.assertEqual(batch.get(t['argv'], '--chr'), t['chromosome'])
            self.assertEqual(batch.get(t['argv'], '--dir-target'), str(self.root / 'gen/1kg/37'))
            script = Path(t['cmd']).read_text()
            self.assertIn('GU_ARG_BATCH_WORKER=1', script)
            self.assertIn('trap finish EXIT', script)
            subprocess.run(['bash', '-n', t['cmd']], check=True)
        self.assertFalse((run / 'submission.jsonl').exists())

    def test_build_submits_twelve(self):
        self.make_plan(submit=True)

    def test_real_submission_journal_with_fake_bsub(self):
        run, plan = self.make_plan()
        responses = [subprocess.CompletedProcess([], 0, f'Job <{200 + i}> is submitted.\n', '')
                     for i in range(12)]
        probe = subprocess.CompletedProcess([], 0, ''.join(f'{200 + i} PEND\n' for i in range(12)), '')
        with patch.dict(os.environ, {'SUBMIT_WAIT': 'FALSE'}), \
             patch.object(batch.submit_jobs, 'bsub_process', side_effect=responses) as bsub, \
             patch.object(batch.subprocess, 'run', return_value=probe):
            batch.manage('submit', run)
            self.assertEqual(bsub.call_count, 12)
            # Repeated submit/retry must never duplicate pending jobs.
            batch.manage('submit', run)
            batch.manage('retry', run)
            self.assertEqual(bsub.call_count, 12)
        journal = [json.loads(x) for x in (run / 'submission.jsonl').read_text().splitlines()]
        accepted = [x for x in journal if x['status'] == 'submitted']
        self.assertEqual(len(accepted), 12)
        self.assertEqual(len({x['job_id'] for x in accepted}), 12)

    def test_partial_submission_does_not_repeat_uncertain_job(self):
        run, plan = self.make_plan()
        responses = [subprocess.CompletedProcess([], 0, 'Job <300> is submitted.\n', ''),
                     subprocess.TimeoutExpired(['bsub'], 1)]
        with patch.dict(os.environ, {'SUBMIT_WAIT': 'FALSE'}), \
             patch.object(batch.submit_jobs, 'bsub_process', side_effect=responses):
            with self.assertRaises(RuntimeError):
                batch.manage('submit', run)
        probe = subprocess.CompletedProcess([], 0, '300 PEND\n', '')
        with patch.object(batch.subprocess, 'run', return_value=probe), \
             patch.object(batch.submit_jobs, 'run') as submit:
            batch.manage('submit', run)
        self.assertEqual(submit.call_args.args[0], [Path(t['cmd']) for t in plan['tasks'][2:]])

    def test_retry_excludes_active_done_unknown_and_unsent(self):
        run, plan = self.make_plan()
        tasks = plan['tasks']
        entries = [dict(cmd=t['cmd'], job_id=str(100 + i), status='submitted')
                   for i, t in enumerate(tasks[:5])]
        # An interrupted/uncertain submission without a known ID must not be retried.
        entries.append(dict(cmd=tasks[5]['cmd'], status='submitting'))
        (run / 'submission.jsonl').write_text(''.join(json.dumps(e) + '\n' for e in entries))
        output = '100 DONE\n101 EXIT\n102 RUN\n103 PEND\n'  # 104 absent => UNKNOWN
        probe = subprocess.CompletedProcess([], 0, stdout=output, stderr='')
        with patch.object(batch.subprocess, 'run', return_value=probe), patch.object(batch.submit_jobs, 'run') as submit:
            batch.manage('retry', run)
            self.assertEqual(submit.call_args.args[0], [Path(tasks[1]['cmd'])])
            submit.reset_mock()
            batch.manage('submit', run)
            self.assertEqual(submit.call_args.args[0], [Path(t['cmd']) for t in tasks[6:]])
            submit.reset_mock()
            batch.manage('retry', run, plan_only=True)
            submit.assert_not_called()
        # Terminal states remain queryable after LSF history expires.
        with patch.object(batch.subprocess, 'run', side_effect=FileNotFoundError):
            rows = batch.task_states(plan)
        self.assertEqual([r['state'] for r in rows[:6]], ['DONE', 'EXIT', 'UNKNOWN', 'UNKNOWN', 'UNKNOWN', 'UNKNOWN'])

    def test_latest_attempt_ignores_old_failure(self):
        run, plan = self.make_plan()
        t = plan['tasks'][0]
        (Path(t['directory']) / 'attempt.100.json').write_text(json.dumps(dict(state='EXIT', elapsed_seconds=20)))
        entries = [dict(cmd=t['cmd'], job_id='100', status='submitted'),
                   dict(cmd=t['cmd'], job_id='101', status='submitted')]
        (run / 'submission.jsonl').write_text(''.join(json.dumps(e) + '\n' for e in entries))
        probe = subprocess.CompletedProcess([], 0, stdout='101 PEND\n', stderr='')
        with patch.object(batch.subprocess, 'run', return_value=probe), patch.object(batch.submit_jobs, 'run') as submit:
            rows = batch.manage('retry', run)
        submit.assert_not_called()
        self.assertEqual(rows[0]['state'], 'PEND')
        self.assertIsNone(rows[0]['elapsed_seconds'])

    def test_timing_and_failure_receipt(self):
        with patch.dict(os.environ, {'LSB_JOBID': '123'}), patch.object(batch.time, 'time', side_effect=[100, 105]):
            batch.record(self.root, 'start')
            batch.record(self.root, 'finish', 7)
        receipt = json.loads((self.root / 'attempt.123.json').read_text())
        self.assertEqual((receipt['state'], receipt['elapsed_seconds'], receipt['exit_code']), ('EXIT', 5, 7))

    def test_invalid_selectors(self):
        for value in ('1,1', '0', '1-23', 'X,23', ''):
            with self.assertRaises(ValueError):
                batch.chromosomes(value)
        self.assertEqual(batch.chromosomes('chr1,20,chrX'), ['1', '20', 'X'])
        with self.assertRaises(ValueError):
            self.make_plan(('--slots', '2'))


if __name__ == '__main__':
    unittest.main()
