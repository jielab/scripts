"""Generate real LD command files without data, LSF, Conda or LD computation."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CommandFiles(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        self.env = dict(os.environ, OUT=str(self.out), PYTHON=sys.executable,
                        CHRS='1,2', GRCH='38', THREADS='16', MEMORY_MB='65536',
                        MEMORY_DEFAULT_GIB='160', LD_MEMORY_GIB='160',
                        LD_WALLTIME='120:00', LD_QUEUE='33', submit='hpc',
                        SUBMIT='hpc', jobs_N='16', SUBMIT_WAIT='FALSE',
                        DIR_WORK='/custom/work', DIR_DAT='/custom/data',
                        DIR_BIG='/custom/big', DIR_TMP='/custom/tmp')
        self.env.pop('LD_SLOTS', None)

    def generate(self):
        return subprocess.run(['bash', str(ROOT / 'f/jobs.sh'), '1'],
                              env=self.env, capture_output=True, text=True)

    def test_compute_and_merge_requests(self):
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        job = (self.out / 'chr1.cmd').read_text()
        for text in ['#BSUB -n 40', '#BSUB -M 160GB', '#BSUB -W 120:00',
                     'rusage[mem=160GB/host]', 'export DIR_WORK=/custom/work',
                     'export DIR_DAT=/custom/data']:
            self.assertIn(text, job)
        merge = (self.out / 'merge.cmd').read_text()
        for text in ['#BSUB -n 10', '#BSUB -M 4GB', '#BSUB -W 1:00']:
            self.assertIn(text, merge)
        for cmd in self.out.glob('*.cmd'):
            self.assertEqual(subprocess.run(['bash', '-n', str(cmd)]).returncode, 0)

    def test_oversized_request_stops_before_submission(self):
        self.env['LD_MEMORY_GIB'] = '896'
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('single-node budget', result.stderr)
        self.assertFalse(list(self.out.glob('*.cmd')))

    def test_auto_queue_and_plink_budget(self):
        self.env.update(LD_QUEUE='auto', LD_MEMORY_GIB='448')
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('#BSUB -q 38', (self.out / 'chr1.cmd').read_text())
        self.env['MEMORY_MB'] = str(512 * 1024)
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('exceeds', result.stderr)

    def test_queue_33_rounds_up_larger_cpu_or_memory_share(self):
        self.env.update(THREADS='4', MEMORY_MB='4096')
        for requested, memory, expected in [(5, 48, 10), (5, 49, 20),
                                             (5, 96, 20), (5, 97, 30),
                                             (5, 144, 30), (5, 145, 40),
                                             (11, 4, 20), (21, 4, 30), (31, 4, 40)]:
            with self.subTest(slots=requested, memory_gib=memory):
                self.env.update(LD_SLOTS=str(requested), LD_MEMORY_GIB=str(memory))
                result = self.generate()
                self.assertEqual(result.returncode, 0, result.stderr)
                job = (self.out / 'chr1.cmd').read_text()
                self.assertIn(f'#BSUB -n {expected}\n', job)
                self.assertIn(f'LSB_DJOB_NUMPROC:-0}} < {expected}', job)

    def test_high_memory_queues_keep_requested_slots(self):
        for queue in ('v4-3t', 'v3-6t'):
            with self.subTest(queue=queue):
                self.env.update(LD_QUEUE=queue, THREADS='4', LD_SLOTS='5')
                result = self.generate()
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('#BSUB -n 5\n', (self.out / 'chr1.cmd').read_text())
                self.assertIn('#BSUB -n 1\n', (self.out / 'merge.cmd').read_text())

    def test_auto_queue_33_applies_rounding_and_small_slot_request_stops(self):
        self.env['LD_QUEUE'] = 'auto'
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        job = (self.out / 'chr1.cmd').read_text()
        self.assertIn('#BSUB -q 33\n', job)
        self.assertIn('#BSUB -n 40\n', job)
        self.env['LD_SLOTS'] = '4'
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('THREADS plus one', result.stderr)

    def test_cli_queue_and_slots_override_environment_and_runtime_defaults(self):
        self.env.update(submit='local', SUBMIT='local', LD_QUEUE='38', LD_SLOTS='20',
                        LD_PYTHON=sys.executable, jobs_N='1')
        (self.out / 'ld.runtime.env').write_text('export LD_QUEUE=v3-6t LD_SLOTS=30\n')
        result = subprocess.run(['bash', str(ROOT / 'ld.sh'), '--grch', '38', '--chr', '1',
                                 '--write-cmd', '-q', '33', '-n', '40'], env=self.env,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        job = (self.out / 'chr1.cmd').read_text()
        self.assertIn('#BSUB -q 33\n', job)
        self.assertIn('#BSUB -n 40\n', job)


if __name__ == '__main__':
    unittest.main()
