"""Exercise real support production, release, reuse and failure boundaries."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from abm_core import pipeline, storage
from abm_core.model import file_sha256


class BusyReasonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.proc = self.root / "321"
        self.proc.mkdir()
        self.stat = self.proc / "stat"
        self.stat.write_text("321 (worker) R")
        (self.proc / "comm").write_text("worker")
        (self.proc / "cmdline").write_bytes(b"python\0worker.py\0")
        self.maps = self.proc / "maps"
        self.maps.write_text("")
        (self.proc / "fd").mkdir()
        self.target = self.root / "geometry.npy"
        self.target.write_bytes(b"persistent cache")
        iterdir = Path.iterdir
        scan = patch.object(Path, "iterdir", lambda path: iterdir(self.root) if path == Path("/proc")
                            else iterdir(path))
        scan.start()
        self.addCleanup(scan.stop)

    def deny_maps(self):
        read = Path.read_text
        def guarded(path, *args, **kwargs):
            if path == self.maps:
                raise PermissionError("process is exiting")
            return read(path, *args, **kwargs)
        return patch.object(Path, "read_text", guarded)

    def test_exit_race_does_not_report_a_cache_user(self):
        for final_state in ("Z", "X", "gone"):
            with self.subTest(state=final_state):
                self.stat.write_text("321 (worker) R")
                def exit_process(_):
                    if final_state == "gone":
                        shutil.rmtree(self.proc)
                    else:
                        self.stat.write_text("321 (worker) " + final_state)
                with self.deny_maps(), patch.object(storage.time, "sleep", exit_process):
                    self.assertIsNone(storage.busy_reason([self.target]))

    def test_disappearing_cmdline_after_permission_denial_is_safe(self):
        (self.proc / "cmdline").unlink()
        with self.deny_maps():
            self.assertIsNone(storage.busy_reason([self.target]))

    def test_slow_exit_race_does_not_report_a_cache_user(self):
        elapsed = 0.
        def exit_process(delay):
            nonlocal elapsed
            elapsed += delay
            if elapsed >= .1:
                self.stat.write_text("321 (worker) Z")
        with self.deny_maps(), patch.object(storage.time, "sleep", exit_process):
            self.assertIsNone(storage.busy_reason([self.target]))
        self.assertGreaterEqual(elapsed, .1)

    def test_exit_race_during_fd_read_does_not_block(self):
        (self.proc / "fd/3").symlink_to(self.target)
        with patch.object(storage.os, "readlink", side_effect=PermissionError("exiting")), \
             patch.object(storage.time, "sleep", lambda _: self.stat.write_text("321 (worker) Z")):
            self.assertIsNone(storage.busy_reason([self.target]))

    def test_persistent_permission_denial_still_blocks(self):
        with self.deny_maps(), patch.object(storage.time, "sleep"):
            self.assertEqual(storage.busy_reason([self.target]), "cannot inspect PID 321")
        self.assertEqual(self.target.read_bytes(), b"persistent cache")

    def test_empty_cache_does_not_require_process_inspection(self):
        self.target.unlink()
        for paths in ([], [self.target], [self.target, self.root / "donor_geometry.npy"]):
            with self.subTest(paths=paths), self.deny_maps(), \
                 patch.object(storage.time, "sleep") as retry:
                self.assertIsNone(storage.busy_reason(iter(paths)))
                retry.assert_not_called()

    def test_missing_scratch_does_not_bypass_existing_donor_protection(self):
        donor = self.root / "donor_geometry.npy"
        self.target.rename(donor)
        with self.deny_maps(), patch.object(storage.time, "sleep"):
            self.assertEqual(storage.busy_reason([self.target, donor]), "cannot inspect PID 321")
        self.assertEqual(donor.read_bytes(), b"persistent cache")

    def test_dangling_symlink_is_not_an_empty_cache(self):
        self.target.unlink()
        self.target.symlink_to(self.root / "missing.npy")
        with self.deny_maps(), patch.object(storage.time, "sleep"):
            self.assertEqual(storage.busy_reason([self.target]), "cannot inspect PID 321")

    def test_retry_still_detects_mapped_and_open_cache(self):
        for kind in ("mapped", "open"):
            with self.subTest(kind=kind):
                denied = self.deny_maps()
                denied.start()
                self.addCleanup(denied.stop)
                def readable(_):
                    denied.stop()
                    if kind == "mapped":
                        self.maps.write_text("0000-1000 rw-s 00000000 00:00 1 " + str(self.target))
                    else:
                        self.maps.write_text("")
                        (self.proc / "fd/3").symlink_to(self.target)
                with patch.object(storage.time, "sleep", readable):
                    reason = storage.busy_reason([self.target])
                self.assertIn(kind, reason)
                self.assertIn("321", reason)


class SupportStorageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cache = Path(self.tmp.name)
        n = 80
        self.frame = pd.DataFrame(dict(id=[f"i{i}" for i in range(n)],
            group_id=[f"g{i}" for i in range(n)],
            role=["build"] * 40 + ["tune_model"] * 20 + ["test"] * 20))
        self.x = np.random.default_rng(9).normal(size=(n, 4, 1))
        self.cfg = json.loads((pipeline.ROOT / "config/abm.json").read_text())
        self.cfg["support"].update(radius_k_ref=5, min_count=1, min_group_ess=1)

    def prepare(self):
        return pipeline.prepare_support(self.frame, self.x, np.arange(80), ["scores"] * 4,
                                        {}, self.cache, self.cfg, "cpu", False)

    def legacy_cache(self):
        with patch.object(storage, "release_verified_scratch", return_value=0):
            self.prepare()
        self.assertTrue((self.cache / "geometry.npy").is_file())

    def test_automatic_release_preserves_frozen_receipt_donor_and_reuse(self):
        first, saved = self.prepare()
        complete_sha = file_sha256(self.cache / "COMPLETE.json")
        release = json.loads((self.cache / storage.RELEASE_FILE).read_text())
        self.assertEqual(release["complete_sha256"], complete_sha)
        self.assertEqual(release["files"]["geometry.npy"]["sha256"], saved["files"]["geometry.npy"])
        self.assertFalse((self.cache / "geometry.npy").exists())
        self.assertFalse((self.cache / "radius_geometry.npy").exists())
        self.assertEqual(file_sha256(self.cache / "donor_geometry.npy"), saved["files"]["donor_geometry.npy"])
        with patch.object(pipeline, "count_twins", side_effect=AssertionError("must reuse")):
            second, reused = self.prepare()
        pd.testing.assert_frame_equal(first, second)
        self.assertTrue(reused["cache_reused"])
        self.assertEqual(file_sha256(self.cache / "COMPLETE.json"), complete_sha)

    def test_legacy_dry_run_apply_and_idempotence(self):
        self.legacy_cache()
        (self.cache / "radius_geometry.npy").write_bytes(b"legacy scratch")
        before = (self.cache / "COMPLETE.json").read_bytes()
        size = storage.clean_cache(self.cache)
        self.assertGreater(size, 0)
        self.assertTrue((self.cache / "geometry.npy").exists())
        self.assertEqual(storage.clean_cache(self.cache, apply=True), size)
        self.assertEqual(storage.clean_cache(self.cache, apply=True), 0)
        self.assertEqual((self.cache / "COMPLETE.json").read_bytes(), before)
        storage.validate_support(self.cache)

    def test_open_or_mapped_files_and_lock_are_retained(self):
        self.legacy_cache()
        with (self.cache / "geometry.npy").open("rb"):
            self.assertEqual(storage.clean_cache(self.cache, apply=True), 0)
        mapped = np.load(self.cache / "geometry.npy", mmap_mode="r")
        try:
            self.assertEqual(storage.clean_cache(self.cache, apply=True), 0)
        finally:
            mapped._mmap.close()
        with storage.support_lock(self.cache):
            with self.assertRaises(BlockingIOError):
                storage.clean_cache(self.cache, apply=True)
        self.assertTrue((self.cache / "geometry.npy").exists())
        self.assertFalse((self.cache / storage.RELEASE_FILE).exists())

    def test_locked_empty_cache_cannot_start_building(self):
        with storage.support_lock(self.cache), patch.object(pipeline, "_prepare_support") as build:
            with self.assertRaises(BlockingIOError):
                self.prepare()
            build.assert_not_called()

    def test_corrupt_retained_output_prevents_deletion(self):
        self.legacy_cache()
        (self.cache / "bank.parquet").write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "Shared support cache changed"):
            storage.clean_cache(self.cache, apply=True)
        self.assertTrue((self.cache / "geometry.npy").exists())

    def test_changed_release_binding_cannot_excuse_missing_scratch(self):
        self.prepare()
        path = self.cache / storage.RELEASE_FILE
        receipt = json.loads(path.read_text())
        receipt["complete_sha256"] = "wrong completion"
        path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "different completion"):
            storage.validate_support(self.cache)

    def test_symlinked_scratch_is_retained(self):
        self.legacy_cache()
        path = self.cache / "geometry.npy"
        target = self.cache / "do-not-delete.npy"
        path.rename(target)
        path.symlink_to(target)
        with self.assertRaisesRegex(ValueError, "Unsafe support cache file"):
            storage.clean_cache(self.cache, apply=True)
        self.assertTrue(path.is_symlink())
        self.assertTrue(target.is_file())

    def test_missing_scratch_requires_matching_release_and_missing_donor_always_fails(self):
        self.legacy_cache()
        (self.cache / "geometry.npy").unlink()
        with self.assertRaisesRegex(ValueError, "Shared support cache changed"):
            self.prepare()
        (self.cache / "COMPLETE.json").unlink()
        self.prepare()
        (self.cache / "donor_geometry.npy").unlink()
        with self.assertRaisesRegex(ValueError, "donor_geometry.npy"):
            self.prepare()

    def test_interrupted_release_retries_without_modifying_completion(self):
        self.legacy_cache()
        unlink = Path.unlink
        def fail_scratch(path, *args, **kwargs):
            if path.name == "geometry.npy":
                raise OSError("injected unlink failure")
            return unlink(path, *args, **kwargs)
        with patch.object(Path, "unlink", fail_scratch):
            with self.assertRaisesRegex(OSError, "injected"):
                storage.clean_cache(self.cache, apply=True)
        storage.validate_support(self.cache)
        self.assertGreater(storage.clean_cache(self.cache, apply=True), 0)

    def test_failures_and_interrupts_close_and_remove_large_incomplete_arrays(self):
        for phase, error in (("estimate_radius", RuntimeError("injected")),
                             ("count_twins", KeyboardInterrupt())):
            with self.subTest(phase=phase), patch.object(pipeline, phase, side_effect=error):
                with self.assertRaises(type(error)):
                    self.prepare()
            self.assertFalse((self.cache / "COMPLETE.json").exists())
            for name in (*storage.SCRATCH, "donor_geometry.npy"):
                self.assertFalse((self.cache / name).exists())

    def test_radius_is_released_before_full_query_count(self):
        count = pipeline.count_twins
        def check(*args, **kwargs):
            self.assertFalse((self.cache / "radius_geometry.npy").exists())
            return count(*args, **kwargs)
        with patch.object(pipeline, "count_twins", side_effect=check):
            self.prepare()

    def test_incomplete_cache_is_not_deleted_by_manual_cleanup(self):
        (self.cache / "geometry.npy").write_bytes(b"unfinished")
        self.assertEqual(storage.clean_cache(self.cache, apply=True), 0)
        self.assertTrue((self.cache / "geometry.npy").exists())

    def analysis_cache(self):
        root = self.cache
        self.cache = root / "cache/abm/support/height/legacy-key/_LOCAL_ONLY"
        self.cache.mkdir(parents=True)
        self.legacy_cache()
        return root

    def test_runner_startup_cleans_legacy_cache_when_phases_are_reused(self):
        root = self.analysis_cache()
        def inputs(*args):
            self.assertFalse((self.cache / "geometry.npy").exists())
        with patch("abm_core.cold_start.ensure_inputs", side_effect=inputs), \
             patch.object(pipeline, "_backend", return_value=("cpu", None)), \
             patch.object(pipeline, "run_trait") as reused:
            pipeline.main(["--analysis-root", str(root), "--traits", "height", "--device", "cpu"])
        reused.assert_called_once()
        storage.validate_support(self.cache)

    def test_runner_exit_retries_busy_cache_without_masking_later_failure(self):
        root = self.analysis_cache()
        with (self.cache / "geometry.npy").open("rb") as reader:
            def later_failure(*args):
                self.assertTrue((self.cache / "geometry.npy").exists())
                reader.close()
                raise RuntimeError("later phase failed")
            with patch("abm_core.cold_start.ensure_inputs"), \
                 patch.object(pipeline, "_backend", return_value=("cpu", None)), \
                 patch.object(pipeline, "run_trait", side_effect=later_failure):
                with self.assertRaisesRegex(RuntimeError, "later phase failed"):
                    pipeline.main(["--analysis-root", str(root), "--traits", "height", "--device", "cpu"])
        self.assertFalse((self.cache / "geometry.npy").exists())
        storage.validate_support(self.cache)

    def test_historical_validation_failure_does_not_block_current_analysis(self):
        root = self.analysis_cache()
        (self.cache / "bank.parquet").write_bytes(b"corrupted historical output")
        with patch("abm_core.cold_start.ensure_inputs"), \
             patch.object(pipeline, "_backend", return_value=("cpu", None)), \
             patch.object(pipeline, "run_trait") as current:
            pipeline.main(["--analysis-root", str(root), "--traits", "height", "--device", "cpu"])
        current.assert_called_once()
        self.assertTrue((self.cache / "geometry.npy").exists())

    def test_sweep_only_cleans_requested_traits(self):
        root = self.analysis_cache()
        self.assertEqual(storage.clean_completed_support(root, traits=["ldl"], apply=True), 0)
        self.assertTrue((self.cache / "geometry.npy").exists())


if __name__ == "__main__":
    unittest.main()
