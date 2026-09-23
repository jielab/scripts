"""Completed runs are kept; failed runs restart while holding an exclusive lock."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "f"))
from common import prepare_run_directory, run_complete, run_lock


class RunLifecycleTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def manifest(self):
        (self.root/"manifest.json").write_text(json.dumps(
            dict(version="5.1.0", config={}, signature="test")))

    def completed(self, train_only=False):
        self.manifest()
        for name in ["MODEL_FROZEN.json", "model_bundle.joblib", "REPORT.md", "test_metrics.csv"]:
            (self.root/name).write_text("test artifact")
        marker = "TRAIN_DONE.json" if train_only else "DONE.json"
        (self.root/marker).write_text('{"version": "5.1.0"}')

    def test_new_directory_starts(self):
        with run_lock(self.root):
            self.assertTrue(prepare_run_directory(self.root))

    def test_complete_run_is_skipped_without_changing_any_artifact(self):
        self.completed()
        before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in self.root.iterdir()}
        with run_lock(self.root):
            self.assertFalse(prepare_run_directory(self.root))
        after = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in self.root.iterdir()}
        self.assertEqual(before, after)

    def test_failed_run_restarts_and_retains_exclusive_lock_during_deletion(self):
        self.manifest()
        (self.root/"input").mkdir()
        (self.root/"input/raw.npy").write_bytes(b"partial")
        with run_lock(self.root):
            self.assertTrue(prepare_run_directory(self.root))
            self.assertEqual([p.name for p in self.root.iterdir()], [".lock"])
            with self.assertRaisesRegex(RuntimeError, "locked"):
                with run_lock(self.root):
                    self.fail("Concurrent training must not acquire the lock")
        self.assertFalse((self.root/".lock").exists())

    def test_frozen_model_without_final_marker_restarts(self):
        self.completed()
        (self.root/"DONE.json").unlink()
        with run_lock(self.root):
            self.assertTrue(prepare_run_directory(self.root))
            self.assertFalse((self.root/"MODEL_FROZEN.json").exists())

    def test_missing_output_or_broken_completion_marker_is_incomplete(self):
        for missing in ["model_bundle.joblib", "REPORT.md", "test_metrics.csv"]:
            with self.subTest(missing=missing):
                self.completed()
                (self.root/missing).unlink()
                self.assertFalse(run_complete(self.root))
        self.completed()
        (self.root/"DONE.json").write_text("{")
        self.assertFalse(run_complete(self.root))

    def test_train_only_requires_training_to_finish_after_freeze(self):
        self.completed(train_only=True)
        self.assertTrue(run_complete(self.root, train_only=True))
        self.assertFalse(run_complete(self.root))
        (self.root/"TRAIN_DONE.json").unlink()
        self.assertFalse(run_complete(self.root, train_only=True))
        self.completed()
        self.assertTrue(run_complete(self.root, train_only=True))

    def test_explicit_replace_rebuilds_completed_run(self):
        self.completed()
        with run_lock(self.root):
            self.assertTrue(prepare_run_directory(self.root, replace=True))
            self.assertFalse((self.root/"DONE.json").exists())

    def test_explicit_resume_preserves_checkpoint(self):
        self.manifest()
        checkpoint = self.root/"checkpoint.pt"
        checkpoint.write_bytes(b"checkpoint")
        with run_lock(self.root):
            self.assertTrue(prepare_run_directory(self.root, resume=True))
        self.assertEqual(checkpoint.read_bytes(), b"checkpoint")

    def test_unmanaged_directory_is_never_cleared(self):
        data = self.root/"unrelated.csv"
        data.write_text("keep me")
        for manifest in [None, "{}", "broken"]:
            if manifest is not None:
                (self.root/"manifest.json").write_text(manifest)
            with run_lock(self.root), self.assertRaisesRegex(ValueError, "manifest"):
                prepare_run_directory(self.root)
            self.assertEqual(data.read_text(), "keep me")

    def test_existing_lock_prevents_replacement(self):
        self.manifest()
        lock = self.root/".lock"
        lock.write_text(str(os.getpid()))
        with self.assertRaisesRegex(RuntimeError, "locked"):
            with run_lock(self.root):
                prepare_run_directory(self.root)
        self.assertTrue((self.root/"manifest.json").is_file())
        self.assertEqual(lock.read_text(), str(os.getpid()))


if __name__ == "__main__":
    unittest.main()
