"""Cross-language path contracts; never read participant data or submit jobs."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]

class PathAdapters(unittest.TestCase):
    def environment(self, mode='local', **extra):
        return dict(PATH=os.environ['PATH'], HOME=os.environ['HOME'], submit=mode, **extra)

    def python(self, code, **env):
        return subprocess.check_output([sys.executable, '-c',
            'import sys; sys.path.insert(0,"0f"); from workflow_paths import configured_path, expand_path_roots; '+code],
            cwd=ROOT, env=self.environment(**env), text=True).strip()

    def test_python_roots_follow_table_and_custom_categories(self):
        result=self.python('print(configured_path("DIR_UKB")); print(configured_path("DIR_GEN", "ukb/37")); print(configured_path("DIR_RUNTIME_TMP"))',
            DIR_DAT='/custom small', DIR_BIG='/custom large', DIR_GEN='/custom genotype', DIR_RUNTIME_TMP='/custom temporary')
        self.assertEqual(result.splitlines(), ['/custom small/phe/ukb','/custom genotype/ukb/37','/custom temporary'])

    def test_python_modes_are_explicit(self):
        self.assertEqual(self.python('print(configured_path("DIR_UKB"))', mode='hpc'), '/work/sph-huangj/phe/ukb')
        self.assertEqual(self.python('print(configured_path("DIR_UKB"))'), '/mnt/d/phe/ukb')

    def test_config_tokens_preserve_explicit_paths(self):
        value=self.python('print(expand_path_roots("${DIR_GEN}/ukb/37")); print(expand_path_roots("/explicit/genotype")); print(expand_path_roots("${HPC_CPU_QUEUE}"))', DIR_GEN='/custom genotype')
        self.assertEqual(value.splitlines(), ['/custom genotype/ukb/37','/explicit/genotype','${HPC_CPU_QUEUE}'])

    def test_grid_checked_in_defaults_expand_without_data_access(self):
        value=self.python('sys.path.insert(0,"grid/f"); from grid_paths import read_config; c=read_config("grid/config/abm.blocks.json"); print(c["score_root"]); print(c["cache_root"])', DIR_BIG='/custom large', DIR_DAT='/custom small')
        self.assertEqual(value.splitlines(), ['/custom large/gen/ukb/pgs','/custom small/analysis/grid/cache/abm/blocks'])

    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_r_reinitialization_tracks_defaults_but_preserves_overrides(self):
        code='''source("0f/0os.R"); os_paths_init("local"); stopifnot(os_path("DIR_UKB")=="/mnt/d/phe/ukb");
Sys.setenv(DIR_DAT="/custom small", DIR_GEN="/custom genotype"); os_paths_init();
stopifnot(os_path("DIR_UKB")=="/custom small/phe/ukb", os_path("DIR_GEN")=="/custom genotype");
Sys.setenv(submit="hpc"); os_paths_init(); stopifnot(os_path("DIR_UKB")=="/work/sph-huangj/phe/ukb",os_path("DIR_GEN")=="/custom genotype"); cat("OK")'''
        value=subprocess.check_output(['Rscript','-e',code],cwd=ROOT,env=self.environment(),text=True)
        self.assertEqual(value,'OK')

if __name__=='__main__': unittest.main()
