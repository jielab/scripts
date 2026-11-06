import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch, MagicMock

spec = importlib.util.spec_from_file_location('submit_jobs', Path(__file__).with_name('submit_jobs.py'))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def script(self, name, text):
        p = self.root / (name + '.cmd')
        p.write_text('#!/bin/bash\nset -euo pipefail\n' + text + '\n')
        return p

    def test_local_pool_cap_and_wait(self):
        payload = self.root / 'work.py'
        payload.write_text('import time,sys\nfrom pathlib import Path\np=Path(sys.argv[1]); a=time.time(); time.sleep(.15); p.write_text(str(a)+" "+str(time.time()))\n')
        paths = [self.script(str(i), f'python3 "{payload}" "{self.root}/{i}.times"') for i in range(6)]
        results = runner.run(paths, 'local', 2)
        self.assertEqual(len(results), 6)
        events = []
        for i in range(6):
            a,b = map(float, (self.root / f'{i}.times').read_text().split())
            events += [(a,1),(b,-1)]
        active = peak = 0
        for _, delta in sorted(events):
            active += delta; peak = max(peak, active)
        self.assertEqual(peak, 2)

    def test_local_failure_stops_new_work_and_merge(self):
        first = self.script('bad', 'exit 7')
        second = self.script('later', f'touch "{self.root}/later"')
        merge = self.script('merge', f'touch "{self.root}/merged"')
        with self.assertRaises(RuntimeError):
            runner.run([first, second], 'local', 1, merge=merge)
        self.assertFalse((self.root / 'later').exists())
        self.assertFalse((self.root / 'merged').exists())
        self.assertIn('"returncode": 7', next(self.root.glob('submission-*.jsonl')).read_text())

    def test_array_index_without_fake_allocation(self):
        p = self.script('array', 'test -z "${LSB_JOBID:-}"; echo "$LSB_JOBINDEX"')
        result = runner.run_local([p], 2, [1,2,3])
        self.assertEqual(len(result), 3)
        self.assertEqual(Path(str(p)+'.3.out').read_text().strip(), '3')

    def test_lsf_concurrency_and_exact_merge_dependencies(self):
        paths = [self.script(str(i), 'true') for i in range(3)]
        merge = self.script('merge', 'true')
        replies = [subprocess.CompletedProcess([], 0, f'Job <{n}> is submitted', '') for n in range(101,105)]
        with patch.object(runner, 'bsub_process', side_effect=replies) as fake:
            self.assertEqual(runner.run(paths, 'hpc', 2, merge=merge), ['101','102','103','104'])
        self.assertEqual(fake.call_args_list[2].args[0], ['bsub','-w','ended(101)'])
        self.assertEqual(fake.call_args_list[3].args[0], ['bsub','-w','done(101) && done(102) && done(103)'])

    def test_disappearance_is_not_success(self):
        with patch.object(runner.subprocess, 'run', return_value=subprocess.CompletedProcess([],1,'','not found')):
            with self.assertRaises(RuntimeError):
                runner.wait_lsf(['123'])

    def test_submission_timeout_default_and_override(self):
        path = self.script('one', 'true')
        for value, expected in ((None, 300), ('600', 600)):
            with self.subTest(value=value), patch.dict(os.environ):
                os.environ.pop('WORKFLOW_BSUB_TIMEOUT_SECONDS', None)
                if value is not None:
                    os.environ['WORKFLOW_BSUB_TIMEOUT_SECONDS'] = value
                child = MagicMock(returncode=0)
                child.communicate.return_value = ('Job <123> is submitted', '')
                with path.open() as handle, patch.object(runner, 'submission_timeout', wraps=runner.submission_timeout) as setting, \
                     patch.object(runner.subprocess, 'Popen') as call:
                    call.return_value.__enter__.return_value = child
                    runner.bsub_process(['bsub'], handle)
                    self.assertEqual(setting(), expected)
                self.assertEqual(call.call_args.kwargs['env']['LC_ALL'], 'C')

    def test_pending_limit_wait_keeps_one_client_and_accepts_job(self):
        path = self.script('one', 'true')
        child = MagicMock(returncode=0)
        pending = b'Group <example>: Pending job threshold reached. Retrying in 60 seconds...\n'
        child.communicate.side_effect = [subprocess.TimeoutExpired(['bsub'], 300, stderr=pending),
                                        ('Job <123> is submitted', pending.decode())]
        with path.open() as handle, patch.object(runner.subprocess, 'Popen') as start, \
             patch.object(runner.time, 'monotonic', side_effect=[0, 0, 301, 301, 302]):
            start.return_value.__enter__.return_value = child
            result = runner.bsub_process(['bsub'], handle, timeout=300)
        self.assertIn('Job <123>', result.stdout)
        start.assert_called_once()
        child.kill.assert_not_called()

    def test_silent_client_is_killed_and_partial_evidence_preserved(self):
        path = self.script('one', 'true')
        child = MagicMock()
        child.communicate.side_effect = [subprocess.TimeoutExpired(['bsub'], 300), ('partial', 'transport stalled')]
        with path.open() as handle, patch.object(runner.subprocess, 'Popen') as start, \
             patch.object(runner.time, 'monotonic', side_effect=[0, 0, 301]):
            start.return_value.__enter__.return_value = child
            with self.assertRaises(subprocess.TimeoutExpired) as error:
                runner.bsub_process(['bsub'], handle, timeout=300)
        child.kill.assert_called_once()
        self.assertEqual(error.exception.stderr, 'transport stalled')

    def test_explicit_post_wait_rejection_refreshes_and_rewinds_script(self):
        path = self.script('one', 'true')
        pending = 'Group <test>: Pending job threshold reached. Retrying in 60 seconds...\n'
        bodies = []
        def client(argv, handle, timeout):
            bodies.append(handle.read())
            if len(bodies) == 1:
                return subprocess.CompletedProcess(argv, 1, '', pending + 'User permission denied. Job not submitted.')
            return subprocess.CompletedProcess(argv, 0, 'Job <123> is submitted', '')
        with path.open() as handle, patch.object(runner, '_bsub_client', side_effect=client):
            result = runner.bsub_process(['bsub'], handle)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(bodies, [path.read_text(), path.read_text()])

    def test_permission_error_or_reported_id_does_not_refresh(self):
        path = self.script('one', 'true')
        pending = 'Group <test>: Pending job threshold reached. Retrying in 60 seconds...\n'
        for response in ('User permission denied. Job not submitted.',
                         pending + 'Job <123> is submitted\nUser permission denied. Job not submitted.'):
            with path.open() as handle, patch.object(runner, '_bsub_client', return_value=
                    subprocess.CompletedProcess(['bsub'], 1, '', response)) as call:
                runner.bsub_process(['bsub'], handle)
            call.assert_called_once()

    def test_submission_timeout_records_uncertainty_without_resubmitting(self):
        paths = [self.script('one', 'true'), self.script('two', 'true')]
        error = subprocess.TimeoutExpired(['bsub'], 300, output=b'Job <123> is submitted', stderr=b'slow response')
        with patch.object(runner, 'bsub_process', side_effect=error) as call:
            with self.assertRaisesRegex(RuntimeError, 'may already be accepted'):
                runner.run(paths, 'hpc')
            self.assertEqual(call.call_count, 1)
        events = [json.loads(s) for s in next(self.root.glob('submission-*.jsonl')).read_text().splitlines()]
        self.assertEqual(events[-1]['status'], 'uncertain')
        self.assertEqual(events[-1]['job_id'], '123')
        self.assertEqual(events[-1]['stderr'], 'slow response')
        self.assertEqual({e['cmd'] for e in events}, {str(paths[0])})

    def test_invalid_timeout_never_submits(self):
        path = self.script('one', 'true')
        for value in ('0', '-1', 'nan', 'inf', 'bad'):
            with self.subTest(value=value), patch.dict(os.environ, {'WORKFLOW_BSUB_TIMEOUT_SECONDS': value}):
                with patch.object(runner, 'bsub_process') as call:
                    with self.assertRaisesRegex(ValueError, 'finite positive'):
                        runner.run([path], 'hpc')
                    call.assert_not_called()
        self.assertFalse(list(self.root.glob('submission-*.jsonl')))

    def test_hpc_more_than_eight_lanes(self):
        paths = [self.script(str(i), 'true') for i in range(22)]
        replies = [subprocess.CompletedProcess([], 0, f'Job <{n}> is submitted', '') for n in range(101, 123)]
        with patch.object(runner, 'bsub_process', side_effect=replies) as fake:
            self.assertEqual(len(runner.run(paths, 'hpc', 12)), 22)
        self.assertEqual(fake.call_args_list[11].args[0], ['bsub'])
        self.assertEqual(fake.call_args_list[12].args[0], ['bsub', '-w', 'ended(101)'])
        with self.assertRaises(ValueError):
            runner.run(paths, 'local', 12)


if __name__ == '__main__':
    unittest.main()
