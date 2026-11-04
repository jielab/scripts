"""Path contract shared by the QiMing entry points (no data access)."""
import os
from pathlib import Path
import subprocess
import unittest

TABLE = Path(__file__).resolve().parents[1] / '0os.sh'

class PathTableTests(unittest.TestCase):
    def shell(self, body, **overrides):
        env = {'PATH': os.environ['PATH'], 'HOME': os.environ['HOME'], **overrides}
        return subprocess.check_output(['bash', '-eu', '-c', 'source "$1"; ' + body, 'paths', str(TABLE)], env=env, text=True).strip()

    def test_qiming_contract(self):
        out = self.shell('printf "%s\\n" "$DIR_WORK" "$DIR_ANNOT" "$DIR_GEN" "$DIR_GWAS" "$DIR_LD" "$DIR_METAG" "$ANALYSIS_ROOT" "$OS_FILES_ROOT"', submit='hpc')
        self.assertEqual(out.splitlines(), ['/work/sph-huangj', '/work/sph-huangj/data/annot', '/work/sph-huangj/data/gen', '/work/sph-huangj/data/gwas', '/work/sph-huangj/data/LD', '/work/sph-huangj/data/metag', '/work/sph-huangj/analysis', '/work/sph-huangj/files'])

    def test_date_does_not_change_persistent_roots(self):
        out = self.shell('DATE_SCRATCH=2026-11-03; source "$1"; printf "%s\\n" "$GRID_OUTPUT_ROOT" "$ALLELE_ARCHAIC_ROOT" "$METAG_STATE_ROOT"', submit='hpc')
        self.assertEqual(out.splitlines(), ['/work/sph-huangj/analysis/grid', '/work/sph-huangj/data/gen', '/work/sph-huangj/analysis/metag'])

    def test_overrides_survive(self):
        out = self.shell('DATE_SCRATCH=2026-11-03; source "$1"; printf "%s:%s" "$DIR_GEN" "$GRID_ANALYSIS_ROOT"', submit='hpc', DIR_GEN='/custom/gen', GRID_ANALYSIS_ROOT='/custom/output')
        self.assertEqual(out, '/custom/gen:/custom/output')

    def test_mode_switch(self):
        out = self.shell('submit=local; source "$1"; printf "%s:%s:%s:%s" "$DIR_WORK" "$DIR_GEN" "$DIR_METAG" "$ANALYSIS_ROOT"', submit='hpc')
        self.assertEqual(out, '/mnt/d:/mnt/f/gen:/mnt/f/metag:/mnt/d/analysis')

    def test_default_hpc_without_scheduler_detection(self):
        out = self.shell('printf "%s:%s:%s:%s" "$submit" "$HPC_GPU_QUEUE" "$HPC_GPU_MEMORY_GIB" "$HPC_WALLTIME"')
        self.assertEqual(out, 'hpc:4a100-40:160:120:00')

    def test_memory_override_and_mode_switch(self):
        out = self.shell('submit=local; source "$1"; printf "%s" "$MEMORY_DEFAULT_GIB"', submit='hpc')
        self.assertEqual(out, '32')
        out = self.shell('submit=local; source "$1"; printf "%s" "$MEMORY_DEFAULT_GIB"', submit='hpc', MEMORY_DEFAULT_GIB='512')
        self.assertEqual(out, '512')

if __name__ == '__main__':
    unittest.main()
