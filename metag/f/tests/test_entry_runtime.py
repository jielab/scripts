"""Regression coverage for execution boundaries; no scheduler jobs are sent."""
import contextlib
import io
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE))
sys.path.insert(0, str(SOURCE.parents[1] / '0f'))
import hpc_common
import hpc_lsf
import submit_jobs


class EntryRuntimeTests(unittest.TestCase):
    def test_install_boolean_matches_other_public_entries(self):
        for value, expected in [('TRUE', 1), ('FALSE', 0)]:
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                cfg = {'work_dir': tmp}
                with patch.object(sys, 'argv', ['install_tools.py', '--submit', value]), \
                     patch.object(hpc_common, 'load_config', return_value=cfg), \
                     patch.object(hpc_lsf, 'write_cmd', return_value=Path(tmp) / 'install.cmd'), \
                     patch.object(hpc_lsf, 'submit_cmd') as submit, \
                     contextlib.redirect_stdout(io.StringIO()):
                    runpy.run_path(str(SOURCE / 'install_tools.py'), run_name='__main__')
                self.assertEqual(submit.call_count, expected)

    def test_install_shell_consumes_shared_options_before_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = {key: str(root / key) for key in hpc_common.PATH_KEYS}
            cfg.update(root=str(root), convert_memory_gib=32)
            config = root / 'config.json'
            config.write_text(json.dumps(cfg))
            result = subprocess.run([
                'bash', str(SOURCE.parent / 'install.sh'), '--submit', 'local',
                '--jobs-N', '2', '--date-scratch', '2026-11-01', '--submit', 'FALSE',
                '--config', str(config)], text=True, capture_output=True,
                env={**os.environ, 'submit': 'local', 'METAG_PYTHON': sys.executable})
            self.assertEqual(result.returncode, 0, result.stderr)
            command = (root / 'work_dir/install.cmd').read_text()
            self.assertIn('export submit=local WORKFLOW_WORKER=1', command)
            self.assertIn('export DATE_SCRATCH=2026-11-01', command)
            self.assertIn(str(config), command)

    def test_shared_execution_code_changes_frozen_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            code = root / 'metag'
            (code / 'f').mkdir(parents=True)
            (code / '0.README.md').write_text('fixture')
            (root / '0f').mkdir()
            for name in ('0os.sh', 'submit_jobs.py'):
                (root / '0f' / name).write_text('first')
            with patch.object(hpc_common, 'CODE_ROOT', code):
                first = hpc_common.code_identity()
                (root / '0f/0os.sh').write_text('second')
                second = hpc_common.code_identity()
                (root / '0f/submit_jobs.py').write_text('second')
                third = hpc_common.code_identity()
            self.assertEqual(len({first, second, third}), 3)

    def _local(self, memory, physical_gib, *, array=None, jobs=1, scope='site_host', cores=1):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            label = 'fixture' + (f'[1-{array}]' if array else '')
            path = root / 'fixture.cmd'
            path.write_text(f'#BSUB -J "{label}"\n#BSUB -n {cores}\n'
                            f'#BSUB -R "span[hosts=1] rusage[mem={memory}]"\nexit 0\n')
            cfg = {'work_dir': str(root), 'lsf_memory_scope': scope}
            def sysconf(name):
                return {'SC_PHYS_PAGES': physical_gib * 2**30 // 4096,
                        'SC_PAGE_SIZE': 4096}[name]
            with patch.object(hpc_lsf.os, 'sysconf', side_effect=sysconf), \
                 patch.dict(os.environ, {'jobs_N': str(jobs)}), \
                 patch.object(submit_jobs, 'run_local', return_value=[{'returncode': 0}]) as run:
                try:
                    return hpc_lsf.local_submit(cfg, path)
                except RuntimeError:
                    run.assert_not_called()
                    raise

    def test_4096_gib_cannot_pass_a_32_gib_local_host(self):
        with self.assertRaisesRegex(RuntimeError, 'exceeds local RAM'):
            self._local('4096GB', 32)

    def test_feasible_local_allocation_is_allowed(self):
        self.assertEqual(self._local('16GB/host', 32)['status'], 'submitted')

    def test_per_task_memory_is_counted_for_all_cpu_slots(self):
        with self.assertRaisesRegex(RuntimeError, 'exceeds local RAM'):
            self._local('8GB', 32, scope='site_task', cores=8)

    def test_array_concurrency_is_included_in_local_memory_budget(self):
        with self.assertRaisesRegex(RuntimeError, 'exceeds local RAM'):
            self._local('16GB', 32, array=4, jobs=4)


if __name__ == '__main__':
    unittest.main()
