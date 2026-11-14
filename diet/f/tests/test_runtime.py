#!/usr/bin/env python3
"""Check local/HPC selection without launching R or submitting any jobs."""
import os
import json
import tempfile
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[2]

class Runtime(unittest.TestCase):
    def probe(self, mode, args=(), overrides=None):
        env = {k: v for k, v in os.environ.items() if not k.startswith(
            ('DIET_', 'LE8_', 'R_', 'PYTHON', 'DIR_', '_OS_DEFAULT_', 'SUBMIT'))
            and k not in ('submit', 'ANALYSIS_ROOT')}
        env.update(submit=mode, SUBMIT=mode, DIET_R_LIBRARY="/nonexistent/diet-test-r-library",
                   DIET_PYTHON_LIBRARY="/nonexistent/diet-test-python-library",
                   DIET_SYSTEM_PREFIX="/nonexistent/diet-test-system-libs",
                   DIET_GRAPHICS_PREFIX="/nonexistent/diet-test-graphics")
        env.update(overrides or {})
        code = '''set -euo pipefail
HELPERS=$1/../0f
ROOT=$1
shift
source "$HELPERS/0os.sh"
workflow_args "$@"
source "$ROOT/f/runtime_env.sh" "${WORKFLOW_ARGS[@]}"
printf '%s\\n' "$submit" "$DIR_FILES/diet" "$ANALYSIS_ROOT/diet" "$DIET_RSCRIPT" "$DIET_PYTHON" "${R_LIBS_USER-unset}"
'''
        return subprocess.check_output(['bash', '-c', code, 'runtime-test', str(ROOT),
             '--submit', mode, *args], env=env, text=True).splitlines()

    def test_local_keeps_system_r_and_user_library(self):
        out = self.probe('local', overrides={'R_LIBS_USER': '/custom/r-library'})
        self.assertEqual(out[:4], ['local', '/mnt/d/files/diet', '/mnt/d/analysis/diet', 'Rscript'])
        self.assertEqual(out[5], '/custom/r-library')

    def test_hpc_uses_os_table_and_isolates_r(self):
        out = self.probe('hpc', overrides={'R_LIBS_USER': '/old/saige'})
        self.assertEqual(out[:3], ['hpc', '/work/sph-huangj/files/diet', '/work/sph-huangj/analysis/diet'])
        self.assertEqual(out[3:5], ['/work/sph-huangj/.conda/envs/le8/bin/Rscript', '/work/sph-huangj/.conda/envs/le8/bin/python3'])
        self.assertEqual(out[5], 'unset')

    def test_private_system_library_is_explicitly_added(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = Path(tmp) / "lib/R/library"
            library.mkdir(parents=True)
            out = self.probe('hpc', overrides={
                'DIET_SYSTEM_PREFIX': tmp, 'R_LIBS_USER': '/old/saige'})
            self.assertEqual(out[5], str(library))

    def test_cli_overrides_runtime_and_root(self):
        out = self.probe('hpc', ['--dir-dat', '/custom/data', '--r-bin=/custom/Rscript', '--python', '/custom/python'],
                         {'DIET_RSCRIPT': '/env/Rscript', 'DIET_PYTHON': '/env/python'})
        self.assertEqual(out[1:5], ['/custom/data/files/diet', '/custom/data/analysis/diet', '/custom/Rscript', '/custom/python'])

    def test_full_dry_run_does_not_shadow_diet_results_module(self):
        out = subprocess.check_output(['bash', str(ROOT / 'diet.sh'), 'all',
            '--cohort', 'ukb', '--study', 'dose', '--core', '--submit', 'local',
            '--python', os.sys.executable, '--dry-run'], text=True)
        self.assertIn('[web/render]', out)
        self.assertIn('[bayes/collect]', out)

    def test_public_entries_write_distinct_lsf_names(self):
        for entry, name in [('1.disco.sh', 'diet.disco'), ('diet.sh', 'diet')]:
            with self.subTest(entry=entry), tempfile.TemporaryDirectory() as tmp:
                env = dict(os.environ, DIET_PYTHON=os.sys.executable, DIET_RSCRIPT='Rscript')
                subprocess.run(['bash', str(ROOT / entry), 'ukb', '--study', 'dose',
                    '--submit', 'hpc', '--dir-output', tmp, '--plan-only'],
                    env=env, text=True, capture_output=True, check=True)
                plan = json.loads(next(Path(tmp).rglob('plan.json')).read_text())
                self.assertEqual(plan['jobs'][0]['job_name'], name)
                self.assertEqual(plan['jobs'][0]['argv'][1], str(ROOT / entry))
                command = next(Path(tmp).rglob('*.cmd')).read_text()
                self.assertIn('#BSUB -J ' + name + '\n', command)
                self.assertFalse(list(Path(tmp).rglob('submission.json')))

    def test_le8_style_runtime_overrides(self):
        out = self.probe('local', overrides={'R_BIN': '/explicit/Rscript', 'PYTHON_BIN': '/explicit/python'})
        self.assertEqual(out[3:5], ['/explicit/Rscript', '/explicit/python'])

    def test_score_and_cohort_lists_survive_scheduler_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, DIET_PYTHON=os.sys.executable, DIET_RSCRIPT='Rscript')
            subprocess.run(['bash', str(ROOT / '1.disco.sh'), 'ukb,nhanes',
                '--score', 'dash,mind,dose', '--replace', 'FALSE', '--submit', 'hpc',
                '--dir-output', tmp, '--plan-only'], env=env, capture_output=True, check=True)
            plan = json.loads(next(Path(tmp).rglob('plan.json')).read_text())
            argv = plan['jobs'][0]['argv']
            self.assertIn('ukb,nhanes', argv)
            self.assertEqual(argv[argv.index('--score') + 1], 'dash,mind,dose')
            self.assertEqual(argv[argv.index('--replace') + 1], 'FALSE')
            self.assertFalse(list(Path(tmp).rglob('submission.json')))

if __name__ == '__main__':
    unittest.main(verbosity=2)
