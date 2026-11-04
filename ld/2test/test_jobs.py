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

    def generate(self):
        return subprocess.run(['bash', str(ROOT / 'f/jobs.sh'), '1'],
                              env=self.env, capture_output=True, text=True)

    def test_compute_and_merge_requests(self):
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        job = (self.out / 'chr1.cmd').read_text()
        for text in ['#BSUB -n 17', '#BSUB -M 160GB', '#BSUB -W 120:00',
                     'rusage[mem=160GB/host]', 'export DIR_WORK=/custom/work',
                     'export DIR_DAT=/custom/data']:
            self.assertIn(text, job)
        merge = (self.out / 'merge.cmd').read_text()
        for text in ['#BSUB -n 1', '#BSUB -M 4GB', '#BSUB -W 1:00']:
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


if __name__ == '__main__':
    unittest.main()
