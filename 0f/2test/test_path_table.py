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
        self.assertEqual(out.splitlines(), ['/work/sph-huangj', '/work/sph-huangj/data/annot', '/work/sph-huangj/data/gen', '/work/sph-huangj/data/gwas', '/work/sph-huangj/data/LD', '/scratch/2026-11-03/sph-huangj/metag', '/scratch/2026-11-03/sph-huangj/analysis', '/work/sph-huangj/data/files'])

    def test_date_propagates(self):
        out = self.shell('DATE_SCRATCH=2026-11-03; source "$1"; printf "%s\\n" "$GRID_OUTPUT_ROOT" "$ALLELE_ARCHAIC_ROOT" "$METAG_STATE_ROOT"', submit='hpc')
        self.assertEqual(out.splitlines(), ['/scratch/2026-11-03/sph-huangj/analysis/grid', '/work/sph-huangj/data/gen', '/scratch/2026-11-03/sph-huangj/metag'])

    def test_overrides_survive(self):
        out = self.shell('DATE_SCRATCH=2026-11-03; source "$1"; printf "%s:%s" "$DIR_GEN" "$GRID_ANALYSIS_ROOT"', submit='hpc', DIR_GEN='/custom/gen', GRID_ANALYSIS_ROOT='/custom/output')
        self.assertEqual(out, '/custom/gen:/custom/output')

    def test_mode_switch(self):
        out = self.shell('submit=local; source "$1"; printf "%s:%s:%s:%s" "$DIR_WORK" "$DIR_GEN" "$DIR_METAG" "$ANALYSIS_ROOT"', submit='hpc')
        self.assertEqual(out, '/mnt/d:/mnt/f/gen:/mnt/f/metag:/mnt/f/analysis')

if __name__ == '__main__':
    unittest.main()
