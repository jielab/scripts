"""Exercise real shell dispatch with stub model work, without datasets or a GPU."""

import hashlib
import os
from pathlib import Path
import selectors
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]


class EntrypointReloadTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="grid-entrypoint-test-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / "f").mkdir()
        for name in ("grid.sh", "grid.upstream.sh", "f/twin_entry.py", "f/prsformer.sh"):
            shutil.copyfile(ROOT / name, self.root / name)
        (self.root / "f/prsformer.py").touch()
        cache = self.root / "cache"
        roster = cache / "roster"
        roster.mkdir(parents=True)
        for name in ("split.rds", "genotypes.npy", "data.tsv.gz", "variants.tsv.gz", "prepare.json"):
            (cache / name).write_text("stub")
        for name in ("common.keep", "families.tsv.gz", "prsformer.split.tsv.gz"):
            (roster / name).write_text("stub")
        self.python = self.root / "stub-python"
        self.python.write_text(f"#!{sys.executable}\n" + '''
import hashlib
import os
import sys

if sys.argv[1] == "-":
    sys.stdin.read()  # Dependency probe or cache lock key; no real imports.
    print(hashlib.sha256(os.environ["GRID_TEST_ROOT"].encode()).hexdigest())
elif sys.argv[2] == "shared-roster":
    print("ROSTER STUB", flush=True)
elif sys.argv[2] == "train":
    print("TRAIN STUB", flush=True)
    if sys.stdin.readline().strip() != "continue":
        raise SystemExit(99)
    raise SystemExit(int(os.environ["GRID_TEST_STATUS"]))
elif sys.argv[2] == "publish":
    print("PUBLISH STUB", flush=True)
''')
        self.python.chmod(0o700)
        key = hashlib.sha256(str(self.root).encode()).hexdigest()
        self.addCleanup(Path(f"/tmp/grid-prsformer-cache-{key}.lock").unlink, missing_ok=True)

    def command(self, standalone):
        common = ["--python", str(self.python), "--cache-dir", str(self.root / "cache")]
        if standalone:
            return ["bash", str(self.root / "f/prsformer.sh"), *common,
                    "--output-root", str(self.root / "out"), "--resume"]
        return ["bash", str(self.root / "grid.sh"), "prsformer", *common,
                "--out-root", str(self.root / "out/abm"),
                "--prsformer-cache-dir", str(self.root / "cache"),
                "--prsformer-python", str(self.python), "--prsformer-arg", "--resume"]

    def run_training(self, standalone=False, edit=False, status=0):
        env = dict(os.environ, SCRIPT_CONSOLE_ACTIVE="1", GRID_PYTHON=sys.executable, GRID_TEST_ROOT=str(self.root),
                   GRID_TEST_STATUS=str(status))
        with subprocess.Popen(self.command(standalone), env=env, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              start_new_session=True) as proc:
            try:
                output = b""
                deadline = time.monotonic() + 10
                with selectors.DefaultSelector() as selector:
                    selector.register(proc.stdout, selectors.EVENT_READ)
                    while b"TRAIN STUB\n" not in output:
                        remaining = deadline - time.monotonic()
                        self.assertTrue(remaining > 0 and selector.select(remaining), output.decode())
                        chunk = os.read(proc.stdout.fileno(), 65536)
                        self.assertTrue(chunk, output.decode())
                        output += chunk
                if edit:
                    script = self.root / ("f/prsformer.sh" if standalone else "grid.upstream.sh")
                    source = script.read_bytes()
                    marker = b"\ncase $MODE in\n" if standalone else b"\nif [[ $MODULE == prsformer ]]; then\n"
                    tail = source[source.rindex(marker) + 1:]
                    # Preserve the inode but grow the file: the old EOF now points
                    # to its final dispatch, reproducing edits during a long run.
                    script.write_bytes(b"#" + b" " * (len(tail) - 2) + b"\n" + source)
                rest, _ = proc.communicate(b"continue\n", timeout=10)
                return proc.returncode, (output + rest).decode()
            finally:
                if proc.poll() is None:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()

    def assert_once(self, standalone, edit, status):
        code, output = self.run_training(standalone, edit, status)
        self.assertEqual(code, status, output)
        self.assertEqual(output.count("TRAIN STUB\n"), 1, output)
        self.assertEqual(output.count("PUBLISH STUB\n"), int(status == 0), output)
        self.assertEqual(output.count("ROSTER STUB\n"), int(not standalone), output)

    def test_shared_training(self):
        self.assert_once(standalone=False, edit=False, status=0)

    def test_shared_training_during_edit(self):
        self.assert_once(standalone=False, edit=True, status=0)

    def test_shared_failure_during_edit(self):
        self.assert_once(standalone=False, edit=True, status=17)

    def test_standalone_training(self):
        self.assert_once(standalone=True, edit=False, status=0)

    def test_standalone_training_during_edit(self):
        self.assert_once(standalone=True, edit=True, status=0)

    def test_standalone_failure_during_edit(self):
        self.assert_once(standalone=True, edit=True, status=17)


if __name__ == "__main__":
    unittest.main()
