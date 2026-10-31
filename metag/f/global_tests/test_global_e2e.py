"""Real public default entry, native build, species output, reuse and cgroup OOM."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import uuid

from test_global_native import _taxonomy, _sequences


@unittest.skipUnless(os.environ.get("METAG_TEST_CGROUP") == "1", "enable real systemd cgroup tests")
class PublicEndToEndTests(unittest.TestCase):
    def test_memory_oom_is_confined_and_child_inherits_limit(self):
        unit = "metag-global-test-" + uuid.uuid4().hex + ".scope"
        sentinel = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
        try:
            code = """import pathlib, subprocess, sys
group = pathlib.Path('/sys/fs/cgroup') / pathlib.Path('/proc/self/cgroup').read_text().strip().split('::/')[1]
assert (group/'memory.max').read_text().strip() == str(128*1024**2)
assert (group/'memory.swap.max').read_text().strip() == '0'
assert (group/'memory.oom.group').read_text().strip() == '1'
print('LIMITS_VERIFIED', flush=True)
subprocess.run([sys.executable, '-c', 'x=bytearray(256*1024**2); print(len(x))'], check=True)
"""
            result = subprocess.run(["systemd-run", "--user", "--scope", "--quiet", "--unit", unit,
                                     "-p", "MemoryMax=128M", "-p", "MemorySwapMax=0",
                                     "-p", "OOMPolicy=kill", sys.executable, "-c", code],
                                    capture_output=True, text=True, timeout=30)
            self.assertIn("LIMITS_VERIFIED", result.stdout)
            self.assertNotEqual(result.returncode, 0)
            self.assertIsNone(sentinel.poll(), "outside process must survive the scoped OOM")
        finally:
            sentinel.terminate()
            sentinel.wait(timeout=10)

    @unittest.skipUnless(os.environ.get("METAG_KRAKEN2"), "requires real pinned Kraken")
    def test_no_module_builds_benchmarks_splits_two_datasets_and_reuses(self):
        with tempfile.TemporaryDirectory(prefix="metag-global-e2e-") as temporary:
            root = Path(temporary)
            refs, templates = _sequences()
            ref, lib, raw = root/"refs", root/"lib", root/"raw"
            ref.mkdir()
            _taxonomy(lib/"metacache/taxonomy")
            inventory = lib/"reference/inventory.tsv"
            inventory.parent.mkdir()
            states, rows = {}, []
            for index, (tid, sequence) in enumerate(refs):
                path = ref/f"ref{index}.fna"
                path.write_text(f">fixture{index}\n{sequence}\n")
                st = path.stat()
                states[path.name] = [st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_dev, st.st_ino]
                rows.append(f"{path.name}\t{path.name}\t{tid}\n")
            inventory.write_text("relative_path\tfilename\ttaxid\n" + "".join(rows))
            st = inventory.stat()
            (inventory.parent/"check-ref.receipt.json").write_text(json.dumps({
                "status": "complete", "references": states, "outputs": {str(inventory):
                [st.st_mode, st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_dev, st.st_ino]}}))
            originals = {}
            for dataset in ("A", "B"):
                (raw/dataset).mkdir(parents=True)
                for mate in (1, 2):
                    payload = bytearray()
                    for index in range(61):
                        sequence = templates[index % 6][mate*3:mate*3+120]
                        payload.extend(f"@read{index}/{mate}\n{sequence}\n+\n{'I'*120}\n".encode())
                    originals[(dataset, mate)] = bytes(payload)
                    with gzip.open(raw/dataset/f"run_{mate}.fastq.gz", "wb") as handle:
                        handle.write(payload)
            config = root/"metag.conf"
            config.write_text("\n".join(f"{key}={value}" for key, value in {
                "DIR_REF": ref, "DIR_LIB": lib, "DIR_RAW": raw, "DIR_OUT": root/"work",
                "DIR_SPLIT": root/"split", "KEEP_REFERENCE_FASTA": "TRUE"}.items()) + "\n")
            command = ["bash", str(Path(__file__).resolve().parents[2]/"3.global.sh"),
                       "--dat-list", "B,A", "--db-gib", "0.001", "--chunk-mib", "1",
                       "--memory-limit-gib", "2", "--memory-reserve-gib", "0",
                       "--min-free-gib", "0", "--threads", "2", "--benchmark-fragments", "17"]
            env = {**os.environ, "METAG_CONFIG": str(config)}
            env.pop("METAG_GLOBAL_MEMORY_SCOPE", None)
            result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Memory scope verified", result.stdout)
            self.assertIn("PASS all global stages", result.stdout)
            self.assertLess(result.stdout.index("[Sample 1/2] B__"), result.stdout.index("[Sample 2/2] A__"))
            benchmark = list((root/"work/runs_global/benchmarks").glob("*/benchmark.json"))
            self.assertEqual(len(benchmark), 1)
            report = json.loads(benchmark[0].read_text())
            self.assertEqual(report["sample"], "B__run__kraken2")
            self.assertEqual(report["query"]["fragments"], 17)
            self.assertIn("-q", report["query"]["command"])
            self.assertEqual(report["query"]["resource_limits"]["memory_max_bytes"], 2*1024**3)
            print("\nE2E default route verified: native masking/build, quick query, "
                  "2 GiB process-tree limit, 2 datasets, lean species output and exact FASTQ preservation",
                  flush=True)
            completed = {}
            for dataset in ("A", "B"):
                target = root/"split"/dataset/"run__kraken2"
                marker = target/"split.done"
                completed[marker] = marker.stat().st_mtime_ns
                saved = json.loads(marker.read_text())
                self.assertEqual(saved["fragments"], 61)
                self.assertFalse(saved["read_reports"])
                self.assertFalse((target/"read-assignments.tsv.gz").exists())
                with (target/"counts.tsv").open() as handle:
                    counts = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(sum(int(row["fragments"]) for row in counts), 61)
                self.assertTrue(any(row["destination"] == "species" for row in counts))
                for mate in (1, 2):
                    expected = originals[(dataset, mate)].splitlines(keepends=True)
                    expected = sorted(b"".join(expected[i:i+4]) for i in range(0, len(expected), 4))
                    observed = []
                    for path in target.rglob(f"PE.R{mate}.fastq.gz"):
                        with gzip.open(path, "rb") as handle:
                            lines = handle.read().splitlines(keepends=True)
                        observed.extend(b"".join(lines[i:i+4]) for i in range(0, len(lines), 4))
                    self.assertEqual(sorted(observed), expected)
            again = subprocess.run(command, env=env, capture_output=True, text=True, timeout=120)
            self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn("REUSE completed global", again.stdout)
            self.assertEqual(again.stdout.count("SKIP completed global sample"), 2)
            for marker, modified in completed.items():
                self.assertEqual(marker.stat().st_mtime_ns, modified)


if __name__ == "__main__":
    unittest.main()
