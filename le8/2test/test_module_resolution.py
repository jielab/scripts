"""Check C1 imports in fresh processes, before abm is cached by other tests."""
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ModuleResolutionTests(unittest.TestCase):
    def test_c1_entrypoints_keep_project_abm(self):
        code = """
from pathlib import Path
import importlib
import sys
root = Path(sys.argv[1])
sys.path.insert(0, str(root / 'f'))
importlib.import_module(sys.argv[2])
import c1_report
import abm
assert Path(abm.__file__).resolve() == root / 'f' / 'abm.py', abm.__file__
assert hasattr(abm, 'Config')
assert callable(abm.uno_c)
"""
        for entrypoint in ("c1_correlate", "c1_evaluation", "hpc_workflow"):
            with self.subTest(entrypoint=entrypoint):
                result = subprocess.run(
                    [sys.executable, "-c", code, str(ROOT), entrypoint],
                    text=True, capture_output=True, timeout=180,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
