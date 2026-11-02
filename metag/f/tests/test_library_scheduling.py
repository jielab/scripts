"""Check user-visible names, queue resources and the conversion dependency."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_cli
import hpc_common
import hpc_reference


class LibrarySchedulingTests(unittest.TestCase):
    def test_prepared_library_uses_short_conversion_and_readable_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = hpc_common.load_config()
            cfg.update(work_dir=str(root), reference_plan=str(root / 'reference.json'))
            tasks = [{'task_id': 1}]
            context = {'tasks': tasks}
            reference = {'status': 'complete', 'identity': hpc_common.fingerprint(context),
                         'context': context, 'tasks': tasks}
            Path(cfg['reference_plan']).write_text(json.dumps(reference))
            with patch.object(hpc_reference, '_load_plan', return_value=reference), \
                 patch.object(hpc_reference, '_prepared_receipt') as verify, \
                 contextlib.redirect_stdout(io.StringIO()):
                path = hpc_cli.generate_library(cfg, prepared=True)
            verify.assert_called_once()
            plan = json.loads(path.read_text())
            self.assertEqual([c['role'] for c in plan['commands']], ['build', 'convert'])
            build, convert = plan['commands']
            self.assertEqual(convert['depends_on'], 'build')
            for command, name, queue, hours in (
                    (build, 'metag_build', 'smp', 168),
                    (convert, 'metag_kunpeng_index', 'short', 72)):
                lines = Path(command['path']).read_text().splitlines()
                self.assertIn(f'#BSUB -J "{name}"', lines)
                self.assertIn(f'#BSUB -q {queue}', lines)
                self.assertIn(f"#BSUB -n {32 if command is build else 40}", lines)
                self.assertIn(f'#BSUB -W {hours}:00', lines)
            self.assertEqual(cfg['build_threads'], 32)
            self.assertLessEqual(cfg['convert_memory_gib'], 168)
