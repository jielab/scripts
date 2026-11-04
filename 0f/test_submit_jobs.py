import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

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
