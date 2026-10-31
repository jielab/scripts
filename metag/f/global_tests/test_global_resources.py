import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import global_resources as resources


class MemoryTests(unittest.TestCase):
    def fixture(self, root):
        proc = root / "proc"
        group = root / "cgroup"
        (proc / "self").mkdir(parents=True)
        group.mkdir()
        (proc / "meminfo").write_text("MemTotal: 104857600 kB\nMemAvailable: 94371840 kB\n")
        (proc / "self/cgroup").write_text("0::/worker\n")
        return proc, group

    def test_cgroup_parent_limits_availability(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc, group = self.fixture(Path(tmp))
            (group / "memory.max").write_text(str(80 * resources.GIB))
            (group / "memory.current").write_text(str(20 * resources.GIB))
            child = group / "worker"
            child.mkdir()
            (child / "memory.max").write_text("max")
            (child / "memory.current").write_text(str(10 * resources.GIB))
            result = resources.memory_status(proc_root=proc, cgroup_root=group)
            self.assertEqual(result["effective_available_bytes"], 60 * resources.GIB)

    def test_host_available_is_a_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc, group = self.fixture(Path(tmp))
            result = resources.memory_status(proc_root=proc, cgroup_root=group)
            self.assertEqual(result["effective_available_bytes"], 90 * resources.GIB)

    def test_insufficient_memory_never_shrinks_budget(self):
        status = {"effective_available_bytes": 75 * resources.GIB}
        with patch.object(resources, "memory_status", return_value=status):
            with self.assertRaisesRegex(RuntimeError, "no automatic budget reduction"):
                resources.require_memory(64 * resources.GIB, 16)

    def test_bad_memory_values_are_rejected(self):
        with self.assertRaises(ValueError):
            resources.require_memory(float("nan"))
        with self.assertRaises(ValueError):
            resources.require_memory(1, -1)


if __name__ == "__main__":
    unittest.main()
