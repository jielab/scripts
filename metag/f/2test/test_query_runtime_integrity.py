"""Regression checks for same-size rewrites around query publication/reuse."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_query
from hpc_common import file_state, sha256


class QueryRuntimeIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        for name in ("output", "work", "reads", "reference", "metadata", "library", "tools"):
            (root / name).mkdir()
        native = root / "tools" / "kraken" / "bin" / "classify"
        native.parent.mkdir(parents=True)
        # A deterministic protocol fixture: exercise the real streaming,
        # compaction, splitting and publication code without compiling Kraken.
        native.write_text("""#!/usr/bin/env python3
import sys
if '-h' in sys.argv:
    print('-H -t -o -p -T -g -Q -S -M')
    sys.exit(0)
while True:
    header = sys.stdin.buffer.readline()
    if not header:
        break
    sequence, plus, quality = [sys.stdin.buffer.readline() for _ in range(3)]
    taxid = 10 if sequence.startswith(b'A') else 20
    print('C\\t%s\\t%d\\t%d\\t%d:1' %
          (header[1:].strip().decode(), taxid, len(sequence.strip()), taxid))
""")
        native.chmod(0o755)
        tax = root / "metadata" / "taxonomy"
        tax.mkdir()
        (tax / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n10\t|\t1\t|\tspecies\t|\n20\t|\t1\t|\tspecies\t|\n")
        (tax / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "10\t|\tA species\t|\t\t|\tscientific name\t|\n"
            "20\t|\tC species\t|\t\t|\tscientific name\t|\n")
        self.db = root / "library" / "db"
        self.db.mkdir()
        for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
            (self.db / name).write_bytes(b"x")
        self.snapshot = {
            "status": "complete", "identity": "protocol-fixture-db",
            "db_dir": str(self.db), "taxonomy_dir": str(tax),
            "files": {p.name: {"state": file_state(p), "sha256": sha256(p)} for p in self.db.iterdir()},
            "taxonomy_files": {p.name: {"state": file_state(p), "sha256": sha256(p)} for p in tax.iterdir()},
        }
        self.cfg = {f"{key}_dir": str(root / key) for key in
                    ("output", "work", "metadata", "library", "tools")}
        self.cfg.update(reference_raw=str(root / "reference"), reads_raw=str(root / "reads"),
                        output_format="directory", min_free_gib=0, classify_threads=1,
                        classifier_engine="kraken", verify_output_hashes=False)
        self.sample = {"sample_id": "fixture/sample", "dataset": "fixture", "sample": "sample",
                       "r1": str(root / "reads" / "sample.fastq"), "r2": None}
        Path(self.sample["r1"]).write_bytes(b"@r1\n" + b"A" * 40 + b"\n+\n" + b"I" * 40 + b"\n")

    def rewrite_same_size(self, path, before, after):
        path = Path(path)
        content = path.read_bytes()
        changed = content.replace(before, after)
        self.assertNotEqual(content, changed)
        self.assertEqual(len(content), len(changed))
        stat = path.stat()
        path.write_bytes(changed)
        # Avoid depending on filesystem timestamp resolution or a test sleep.
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000000000))

    def run_sample(self):
        return hpc_query.run_sample(self.cfg, self.sample, self.snapshot)

    def test_same_size_raw_replacement_never_skips_old_result(self):
        completed = self.run_sample()
        counts = Path(completed["output"]) / "counts.tsv"
        original_counts = counts.read_bytes()
        old_size_state = file_state(self.sample["r1"])
        self.rewrite_same_size(self.sample["r1"], b"A" * 40, b"C" * 40)
        self.assertEqual(file_state(self.sample["r1"]), old_size_state)
        with self.assertRaisesRegex(RuntimeError, "Output completion identity conflicts"):
            self.run_sample()
        self.assertEqual(counts.read_bytes(), original_counts)

    def test_same_size_published_corruption_fails_without_full_hash_option(self):
        completed = self.run_sample()
        counts = Path(completed["output"]) / "counts.tsv"
        self.rewrite_same_size(counts, b"A species", b"Z species")
        with self.assertRaisesRegex(RuntimeError, "Completed output changed"):
            self.run_sample()

    def test_legacy_output_without_runtime_guards_requires_content_hash(self):
        completed = self.run_sample()
        output = Path(completed["output"])
        marker = output / "complete.json"
        receipt = json.loads(marker.read_text())
        for item in receipt["outputs"].values():
            item.pop("runtime_state")
        marker.write_text(json.dumps(receipt))
        self.assertIsNotNone(hpc_query._complete(output, completed["identity"], verify_hashes=False))
        self.rewrite_same_size(output / "counts.tsv", b"A species", b"Z species")
        with self.assertRaisesRegex(RuntimeError, "Completed output hash mismatch"):
            hpc_query._complete(output, completed["identity"], verify_hashes=False)

    def test_same_size_input_rewrite_during_query_stops_publication(self):
        query_once = hpc_query.query_once

        def query_then_rewrite(*args, **kwargs):
            query = query_once(*args, **kwargs)
            self.rewrite_same_size(self.sample["r1"], b"A" * 40, b"C" * 40)
            return query

        with patch.object(hpc_query, "query_once", side_effect=query_then_rewrite):
            with self.assertRaisesRegex(RuntimeError, "Original FASTQ runtime state changed"):
                self.run_sample()
        self.assertFalse((Path(self.cfg["output_dir"]) / "fixture" / "sample").exists())

    def test_same_size_index_rewrite_during_query_stops_publication(self):
        query_once = hpc_query.query_once

        def query_then_rewrite(*args, **kwargs):
            query = query_once(*args, **kwargs)
            self.rewrite_same_size(self.db / "hash.k2d", b"x", b"y")
            return query

        with patch.object(hpc_query, "query_once", side_effect=query_then_rewrite):
            with self.assertRaisesRegex(RuntimeError, "Database or taxonomy runtime state changed"):
                self.run_sample()
        self.assertFalse((Path(self.cfg["output_dir"]) / "fixture" / "sample").exists())


if __name__ == "__main__":
    unittest.main()
