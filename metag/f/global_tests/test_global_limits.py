import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import global_limits as limits


class ScopeTests(unittest.TestCase):
    def setUp(self):
        self.args = SimpleNamespace(command="all", memory_limit_gib=40,
                                    db_gib=24, memory_reserve_gib=8)
        self.enterContext(patch.dict(os.environ, {limits.MARKER: ""}))

    def test_budget_cannot_exceed_scope(self):
        self.args.db_gib = 32
        with patch.object(limits.subprocess, "Popen") as spawn:
            with self.assertRaisesRegex(RuntimeError, "must exceed"):
                limits.ensure_memory_scope(self.args, [])
        spawn.assert_not_called()

    def test_low_host_headroom_refuses_before_launch(self):
        with patch.object(limits, "memory_status", return_value={"effective_available_bytes": 47*limits.GIB}), \
                patch.object(limits.subprocess, "Popen") as spawn:
            with self.assertRaisesRegex(RuntimeError, "Other tasks are left running"):
                limits.ensure_memory_scope(self.args, [])
        spawn.assert_not_called()

    def test_missing_systemd_never_falls_back_to_unlimited(self):
        with patch.object(limits, "memory_status", return_value={"effective_available_bytes": 60*limits.GIB}), \
                patch.object(limits.shutil, "which", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "refusing"):
                limits.ensure_memory_scope(self.args, [])

    def test_watchdog_stops_only_owned_scope(self):
        child = Mock()
        child.poll.return_value = None
        with patch.object(limits, "memory_status", side_effect=[
                {"effective_available_bytes": 60*limits.GIB},
                {"effective_available_bytes": 7*limits.GIB}]), \
                patch.object(limits.shutil, "which", return_value="/usr/bin/systemd-run"), \
                patch.object(limits.subprocess, "Popen", return_value=child) as spawn, \
                patch.object(limits, "stop_scope") as stop:
            with self.assertRaisesRegex(RuntimeError, "stopping only this"):
                limits.ensure_memory_scope(self.args, ["--dat-list", "A,B"])
        command = spawn.call_args.args[0]
        unit = command[command.index("--unit") + 1]
        stop.assert_called_once_with(unit, child)
        self.assertTrue(unit.startswith("metag-global-"))
        self.assertIn("MemoryMax=" + str(40*limits.GIB), command)
        self.assertIn("MemorySwapMax=0", command)
        self.assertEqual(command[-2:], ["--dat-list", "A,B"])

    def test_forged_guard_and_missing_kernel_limit_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(limits, "current_group", return_value=root):
                with self.assertRaisesRegex(RuntimeError, "does not match"):
                    limits.verify_scope("metag-global-forged.scope", 123)
            unit = "metag-global-test.scope"
            group = root/unit
            group.mkdir()
            (group/"memory.max").write_text("123")
            (group/"memory.swap.max").write_text("0")
            (group/"memory.oom.group").write_text("0")
            with patch.object(limits, "current_group", return_value=group):
                with self.assertRaisesRegex(RuntimeError, "were not applied"):
                    limits.verify_scope(unit, 123)


if __name__ == "__main__":
    unittest.main()
