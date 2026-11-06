"""Exercise bundled public workers, failure propagation and shared preparation."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

GRID = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(GRID / "f"))


class TraitWorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="grid trait worker ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = dict(os.environ, submit="local", SUBMIT="local", WORKFLOW_WORKER="1",
                        WORKFLOW_RESOURCE_MANAGER="local", WORKFLOW_PYTHON=sys.executable,
                        DIR_WORK=str(self.root), DIR_DATA=str(self.root / "data"))

    def step(self, number, code=0):
        program = """import json,os,sys
from pathlib import Path
assert 'GRID_TRAIT_STEPS' not in os.environ
assert os.environ['WORKFLOW_WORKER']=='1'
p=Path(sys.argv[1]); n=int(sys.argv[2]); code=int(sys.argv[3])
seen=json.loads(p.read_text()) if p.exists() else []
assert len(seen)==n-1
seen.append({'step':n,'threads':os.environ['GRID_THREADS']})
p.write_text(json.dumps(seen))
sys.exit(code)
"""
        return dict(key=f"height.step{number}", argv=[sys.executable, "-c", program,
                    str(self.root / "events.json"), str(number), str(code)],
                    env={"GRID_THREADS": str(number)})

    def run_entry(self, filename, steps):
        env = dict(self.env, GRID_TRAIT_STEPS=json.dumps(steps))
        return subprocess.run(["bash", str(GRID/filename), "all", "--traits", "height"],
                              env=env, text=True, capture_output=True, timeout=20)

    def test_all_three_public_entries_run_stages_once_in_order(self):
        for filename in ("2.csx.sh", "3.disco.sh", "4.prsformer.sh"):
            with self.subTest(entry=filename):
                (self.root/"events.json").unlink(missing_ok=True)
                result = self.run_entry(filename, [self.step(1), self.step(2), self.step(3)])
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                events = json.loads((self.root/"events.json").read_text())
                self.assertEqual(events, [dict(step=i, threads=str(i)) for i in (1, 2, 3)])
                self.assertIn("[trait] DONE 3/3", result.stdout)
                self.assertNotIn("[PLAN]", result.stdout)

    def test_failed_stage_stops_the_entire_trait_job(self):
        result = self.run_entry("2.csx.sh", [self.step(1), self.step(2, 7), self.step(3)])
        self.assertEqual(result.returncode, 7, result.stdout+result.stderr)
        self.assertEqual(len(json.loads((self.root/"events.json").read_text())), 2)
        self.assertIn("later stages skipped", result.stderr)
        self.assertNotIn("START 3/3", result.stdout)

    def test_direct_hpc_bundle_requires_a_real_allocation(self):
        env = {k:v for k,v in self.env.items() if not k.startswith("LSB_")}
        env.update(submit="hpc", GRID_TRAIT_STEPS=json.dumps([self.step(1)]))
        result = subprocess.run([sys.executable, str(GRID/"f/trait_workflow.py")],
                                env=env, text=True, capture_output=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("allocated LSF job", result.stderr)
        self.assertFalse((self.root/"events.json").exists())

    def test_concurrent_shared_preparation_waits_and_reuses_the_cache(self):
        from shared_cohort import ensure_inputs
        cache = self.root / "shared-cache"
        cache.mkdir()
        options = SimpleNamespace(cache_dir=str(cache))
        args = SimpleNamespace(analysis_root=str(self.root), complete_features=False)
        started = threading.Event()
        events = []

        def prepare(_):
            receipt = cache / "prepare.json"
            if receipt.exists():
                events.append("reuse")
                return
            events.append("produce")
            started.set()
            time.sleep(.2)
            receipt.write_text('{"options": {}}')

        # Retain ensure_inputs' real run_lock and analysis_guard. Replace only
        # genotype processing and its additional independent file locks.
        with patch("shared_cohort.preparation_command", return_value=["prepare"]), \
             patch("shared_cohort.preparation_options", return_value=(options, {})), \
             patch("shared_cohort.prepare_frozen_roster"), \
             patch("prsformer.prepare", side_effect=prepare), \
             patch("dosage_storage.dosage_lock", side_effect=lambda *a, **kw: nullcontext()), \
             patch("dosage_storage.prsformer_cache_lock", side_effect=lambda *a, **kw: nullcontext()), \
             ThreadPoolExecutor(max_workers=2) as pool:
            first = pool.submit(ensure_inputs, args, {"test": True})
            self.assertTrue(started.wait(10))
            second = pool.submit(ensure_inputs, args, {"test": True})
            first.result(timeout=10)
            second.result(timeout=10)
        self.assertEqual(events, ["produce", "reuse"])
        receipt = json.loads((self.root/"cache/shared_preparation/PREPARED.json").read_text())
        self.assertEqual(receipt["status"], "completed")

    def test_verified_cache_reuse_does_not_request_a_writer_lock(self):
        from shared_cohort import ensure_inputs
        cache = self.root / "shared-cache"
        cache.mkdir()
        (cache / "prepare.json").write_text('{}')
        (cache / "genotypes.npy").write_bytes(b"fixture")
        options = SimpleNamespace(cache_dir=str(cache), dir_gen="fixture", chrs="1", replace=False)
        args = SimpleNamespace(analysis_root=str(self.root), complete_features=False)
        locks = []

        def lock(*args, **kwargs):
            locks.append(kwargs.get("exclusive", False))
            return nullcontext()

        with patch("shared_cohort.preparation_command", return_value=["prepare"]), \
             patch("shared_cohort.preparation_options", return_value=(options, {})), \
             patch("shared_cohort.prepare_frozen_roster"), \
             patch("prsformer.genotype_files", return_value=[]), \
             patch("prsformer.identity", return_value=("verified", {})), \
             patch("prsformer.validate_cache", return_value=True) as validate, \
             patch("prsformer.prepare") as prepare, \
             patch("dosage_storage.dosage_lock", side_effect=lock), \
             patch("dosage_storage.prsformer_cache_lock", side_effect=lock):
            ensure_inputs(args, {"test": True})
        validate.assert_called_once_with(cache, "verified")
        prepare.assert_not_called()
        self.assertEqual(locks, [False, False])


if __name__ == "__main__":
    unittest.main()
