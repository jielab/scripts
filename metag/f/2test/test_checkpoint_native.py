"""Real native interruption/integrity tests (no HPC throughput claims).

METAG_NATIVE_BUILD=/absolute/path/metag_build python -m unittest discover \
  -s metag/f/tests -p 'test_checkpoint_native.py' -v

The tests generate synthetic FASTA/taxonomy in TemporaryDirectory, invoke the
actual C++ executable, deliver real signals, and simulate a partial disk write
with RLIMIT_FSIZE. No production references/checkpoints are touched.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import random
import resource
import shutil
import signal
import subprocess
import tempfile
import unittest

BINARY = os.environ.get("METAG_NATIVE_BUILD")
OFFICIAL = os.environ.get("METAG_OFFICIAL_BUILD")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@unittest.skipUnless(BINARY, "Set METAG_NATIVE_BUILD to the compiled native executable")
class NativeCheckpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = Path(BINARY).resolve()
        cls.fixture = tempfile.TemporaryDirectory(prefix="metag-native-inputs-")
        cls.data = Path(cls.fixture.name)
        (cls.data / "tax").mkdir()
        nodes = [(1, 1, "no rank"), (2, 1, "superkingdom"), (3, 2, "genus"),
                 (10, 3, "species"), (20, 3, "species")]
        (cls.data / "tax/nodes.dmp").write_text("".join(
            f"{t}\t|\t{p}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            for t, p, rank in nodes))
        (cls.data / "tax/names.dmp").write_text("".join(
            f"{t}\t|\tnode{t}\t|\t\t|\tscientific name\t|\n" for t, _, _ in nodes))
        cls.paths = []
        rng = random.Random(412)
        for i, length in enumerate((120000, 1500000, 120000)):
            path = cls.data / f"part{i}.fna"
            seq = "".join(rng.choices("ACGT", k=length))
            path.write_text(f">ref{i}\n" + "\n".join(seq[j:j+80] for j in range(0, length, 80)) + "\n")
            cls.paths.append((path, length))
        (cls.data / "map.tsv").write_text("ref0\t10\nref1\t20\nref2\t10\n")
        cls.manifest = cls.data / "inputs.tsv"
        cls.manifest.write_text("".join(
            f"{p}\t{p.stat().st_size}\t{digest(p)}\t1\t{n}\n" for p, n in cls.paths))
        cls.context = "a" * 64
        cls.capacity = 1200000
        cls.baseline = cls.data / "baseline"
        cls.baseline.mkdir()
        result = subprocess.run(cls.command(cls.baseline), text=True, capture_output=True, env=cls.environment())
        if result.returncode:
            raise RuntimeError(result.stderr)
        cls.expected = digest(cls.baseline / "hash.k2d")

    @classmethod
    def tearDownClass(cls):
        cls.fixture.cleanup()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="metag-native-checkpoint-")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    @classmethod
    def environment(cls):
        return dict(os.environ, OMP_NUM_THREADS="2", K2_DB_READ_THREADS="2")

    @classmethod
    def command(cls, root, *extra):
        return [str(cls.binary), "-H", str(root / "hash.k2d"), "-t", str(root / "taxo.k2d"),
                "-o", str(root / "opts.k2d"), "-m", str(cls.data / "map.tsv"),
                "-n", str(cls.data / "tax"), "-k", "35", "-l", "31", "-p", "2",
                "-B", "16384", "-b", "512", "-c", str(cls.capacity),
                "--metag-inputs", str(cls.manifest), "--metag-context", cls.context,
                "--metag-checkpoint-dir", str(root / "checkpoints"), *map(str, extra)]

    def run_native(self, *extra, **kwargs):
        return subprocess.run(self.command(self.root, *extra), text=True, capture_output=True,
                              env=self.environment(), timeout=90, **kwargs)

    def receipt(self):
        name = (self.root / "checkpoints/CURRENT").read_text().split("\t")[0]
        return json.loads((self.root / "checkpoints" / name / "receipt.json").read_text())

    def pause_first(self):
        result = self.run_native("--metag-stop-after-files", 1)
        self.assertEqual(result.returncode, 75, result.stderr)
        self.assertEqual(self.receipt()["completed_files"], 1)
        self.assertFalse((self.root / "hash.k2d").exists())
        return (self.root / "checkpoints/CURRENT").read_bytes()

    def test_resume_is_byte_identical_skips_completed_input_and_hardlinks_final(self):
        self.pause_first()
        result = self.run_native()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("METAG_RESUME files=1/3", result.stderr)
        self.assertNotIn("METAG_FILE_BEGIN 1/3", result.stderr)
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)
        receipt = self.receipt()
        self.assertTrue(receipt["complete"])
        self.assertEqual(receipt["sequences"], 3)
        self.assertEqual(receipt["bases"], sum(n for _, n in self.paths))
        current_hash = self.root / "checkpoints" / receipt["generation"] / "hash.k2d"
        self.assertEqual(current_hash.stat().st_ino, (self.root / "hash.k2d").stat().st_ino)
        self.assertEqual(len(list((self.root / "checkpoints").glob("cp-*"))), 1)
        self.assertEqual((self.root / "opts.k2d").read_bytes(), (self.baseline / "opts.k2d").read_bytes())
        again = self.run_native()
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertIn("METAG_COMPLETE_REUSE", again.stderr)

    def test_short_write_preserves_previous_current_and_resumes(self):
        previous = self.pause_first()
        def limit_file_size():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (65536, 65536))
        failed = self.run_native(preexec_fn=limit_file_size)
        self.assertNotEqual(failed.returncode, 0)
        self.assertNotEqual(failed.returncode, 75)
        self.assertEqual((self.root / "checkpoints/CURRENT").read_bytes(), previous)
        stages = list((self.root / "checkpoints").glob("stage-*/hash.k2d"))
        self.assertEqual(len(stages), 1)
        self.assertEqual(stages[0].stat().st_size, 65536)
        resumed = self.run_native()
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)
        self.assertFalse(list((self.root / "checkpoints").glob("stage-*")))

    def test_real_sigkill_during_next_input_keeps_committed_progress(self):
        previous = self.pause_first()
        proc = subprocess.Popen(self.command(self.root), text=True, stdout=subprocess.DEVNULL,
                                stderr=subprocess.PIPE, env=self.environment())
        log = []
        try:
            for line in proc.stderr:
                log.append(line)
                if "METAG_FILE_BEGIN 2/3" in line:
                    proc.kill()
                    break
            proc.communicate(timeout=30)
            self.assertEqual(proc.returncode, -signal.SIGKILL, "".join(log))
        finally:
            if proc.poll() is None:
                proc.kill(); proc.wait()
        self.assertEqual((self.root / "checkpoints/CURRENT").read_bytes(), previous)
        resumed = self.run_native()
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertNotIn("METAG_FILE_BEGIN 1/3", resumed.stderr)
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)

    def test_real_sigusr2_finishes_current_file_then_durable_pause(self):
        self.pause_first()
        proc = subprocess.Popen(self.command(self.root), text=True, stdout=subprocess.DEVNULL,
                                stderr=subprocess.PIPE, env=self.environment())
        log = []
        try:
            for line in proc.stderr:
                log.append(line)
                if "METAG_FILE_BEGIN 2/3" in line:
                    proc.send_signal(signal.SIGUSR2)
                    break
            _, remainder = proc.communicate(timeout=60)
            log.append(remainder)
            self.assertEqual(proc.returncode, 75, "".join(log))
        finally:
            if proc.poll() is None:
                proc.kill(); proc.wait()
        self.assertEqual(self.receipt()["completed_files"], 2)
        self.assertEqual(self.receipt()["signal"], signal.SIGUSR2)
        resumed = self.run_native()
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)

    def test_tampered_payload_fails_closed(self):
        previous = self.pause_first()
        receipt = self.receipt()
        path = self.root / "checkpoints" / receipt["generation"] / "hash.k2d"
        with path.open("r+b") as handle:
            handle.seek(-1, 2); byte = handle.read(1)
            handle.seek(-1, 2); handle.write(bytes([byte[0] ^ 1]))
        rejected = self.run_native()
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("checkpoint hash SHA256 mismatch", rejected.stderr)
        self.assertEqual((self.root / "checkpoints/CURRENT").read_bytes(), previous)

    def test_semantic_change_and_fast_mode_rejected(self):
        previous = self.pause_first()
        for extra in (("-p", "1"), ("-b", "256"), ("-F",)):
            rejected = self.run_native(*extra)
            self.assertNotEqual(rejected.returncode, 0, rejected.stderr)
            self.assertEqual((self.root / "checkpoints/CURRENT").read_bytes(), previous)

    def test_recover_first_generation_pointer_publication_gap(self):
        self.pause_first()
        (self.root / "checkpoints/CURRENT").unlink()
        resumed = self.run_native()
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertIn("METAG_RECOVER_FIRST_CURRENT", resumed.stderr)
        self.assertNotIn("METAG_FILE_BEGIN 1/3", resumed.stderr)
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)

    @unittest.skipUnless(OFFICIAL, "Set METAG_OFFICIAL_BUILD for upstream deterministic comparison")
    def test_same_hash_as_original_uninterrupted_deterministic_builder(self):
        command = self.command(self.root)
        end = command.index("--metag-inputs")
        command = [str(Path(OFFICIAL).resolve())] + command[1:end]
        fasta = b"".join(path.read_bytes() for path, _ in self.paths)
        result = subprocess.run(command, input=fasta, capture_output=True,
                                env=self.environment(), timeout=90)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)

    def test_independent_category_builds_keep_identical_taxonomy_and_options(self):
        outputs = []
        for i in (0, 2):
            folder = self.root / f"category{i}"
            folder.mkdir()
            path, bases = self.paths[i]
            manifest = folder / "category.tsv"
            manifest.write_text(f"{path}\t{path.stat().st_size}\t{digest(path)}\t1\t{bases}\n")
            command = self.command(folder)
            command[command.index("--metag-inputs") + 1] = str(manifest)
            command[command.index("-c") + 1] = str(100000 + i*10000)
            result = subprocess.run(command, capture_output=True, text=True,
                                    env=self.environment(), timeout=90)
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(folder)
        for name in ("taxo.k2d", "opts.k2d"):
            self.assertEqual((outputs[0] / name).read_bytes(), (outputs[1] / name).read_bytes())

    def test_changed_completed_input_is_detected_on_resume(self):
        copied = self.root / "first.fna"
        copied.write_bytes(self.paths[0][0].read_bytes())
        manifest = self.root / "inputs.tsv"
        manifest.write_text(self.manifest.read_text().replace(str(self.paths[0][0]), str(copied)))
        command = self.command(self.root)
        command[command.index("--metag-inputs") + 1] = str(manifest)
        first = subprocess.run(command + ["--metag-stop-after-files", "1"], capture_output=True,
                               text=True, env=self.environment(), timeout=90)
        self.assertEqual(first.returncode, 75, first.stderr)
        previous = (self.root / "checkpoints/CURRENT").read_bytes()
        content = copied.read_text()
        offset = content.index("\n") + 1
        content = content[:offset] + ("A" if content[offset] != "A" else "C") + content[offset+1:]
        copied.write_text(content)
        rejected = subprocess.run(command, capture_output=True, text=True,
                                  env=self.environment(), timeout=90)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("input SHA256 mismatch", rejected.stderr)
        self.assertEqual((self.root / "checkpoints/CURRENT").read_bytes(), previous)

    def test_hash_checker_works_when_fork_and_handwritten_vfork_are_forbidden(self):
        """Simulate a huge parent's fork failure; require safe libc spawn instead.

        The shim rejects any direct fork/vfork call and checks the GNU spawn
        flag needed to avoid fork on older glibc. The real full hash build and
        checkpoint resume must still match the deterministic baseline.
        """
        compiler = shutil.which("gcc") or shutil.which("cc")
        if not compiler:
            self.skipTest("C compiler required for fork-failure preload regression")
        shim = self.root / "spawn_probe.c"
        shim.write_text(r'''
#define _GNU_SOURCE
#include <spawn.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <dlfcn.h>

pid_t fork(void) { errno = ENOMEM; return -1; }
pid_t vfork(void) { errno = ENOMEM; return -1; }

int posix_spawnp(pid_t *pid, const char *file,
                 const posix_spawn_file_actions_t *actions,
                 const posix_spawnattr_t *attributes,
                 char *const argv[], char *const envp[]) {
  typedef int (*spawn_fn)(pid_t *, const char *, const posix_spawn_file_actions_t *,
                         const posix_spawnattr_t *, char *const[], char *const[]);
  spawn_fn actual = (spawn_fn) dlsym(RTLD_NEXT, "posix_spawnp");
  if (!actual) return ENOSYS;
  short flags = 0;
  if (!attributes || posix_spawnattr_getflags(attributes, &flags)) return EINVAL;
#ifdef POSIX_SPAWN_USEVFORK
  if (!(flags & POSIX_SPAWN_USEVFORK)) return ENOMEM;
#endif
  const char *marker = getenv("METAG_TEST_SPAWN_MARKER");
  if (marker) {
    int fd = open(marker, O_CREAT | O_WRONLY | O_APPEND, 0600);
    if (fd < 0) return errno;
    const char line[] = "safe-spawn\n";
    if (write(fd, line, sizeof(line)-1) != sizeof(line)-1) { close(fd); return EIO; }
    close(fd);
  }
  return actual(pid, file, actions, attributes, argv, envp);
}
''')
        library = self.root / "spawn_probe.so"
        compiled = subprocess.run([compiler, "-shared", "-fPIC", str(shim), "-ldl", "-o", str(library)],
                                  text=True, capture_output=True, timeout=30)
        self.assertEqual(compiled.returncode, 0, compiled.stderr)
        marker = self.root / "spawns.log"
        environment = self.environment()
        environment["LD_PRELOAD"] = str(library) + ((":" + environment["LD_PRELOAD"]) if environment.get("LD_PRELOAD") else "")
        environment["METAG_TEST_SPAWN_MARKER"] = str(marker)
        first = subprocess.run(self.command(self.root, "--metag-stop-after-files", 1),
                               env=environment, text=True, capture_output=True, timeout=90)
        self.assertEqual(first.returncode, 75, first.stderr)
        second = subprocess.run(self.command(self.root), env=environment, text=True,
                                capture_output=True, timeout=90)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(digest(self.root / "hash.k2d"), self.expected)
        calls = marker.read_text().splitlines()
        self.assertGreater(len(calls), 10)
        self.assertEqual(set(calls), {"safe-spawn"})


if __name__ == "__main__":
    unittest.main()
