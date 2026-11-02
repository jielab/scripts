"""Integration tests for real preparation, resumable category builds and publish.

Set METAG_TEST_TOOLS to an existing synthetic/developer tools directory with the
pinned official Kraken installation and the current private checkpoint builder.
All references, taxonomy, prepared files and indexes are created in temporary
directories. No cluster, production data, or external network is used.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_reference
import partitioned_build as build
import partitioned_reference as prepare
from hpc_common import sha256

TOOLS = os.environ.get("METAG_TEST_TOOLS")


@unittest.skipUnless(TOOLS, "Set METAG_TEST_TOOLS for real category-build integration")
class PartitionedBuildNativeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="metag-category-build-")
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.cfg = {key: str(self.base / name) for key, name in (
            ("reference_raw", "raw"), ("metadata_dir", "metadata"),
            ("library_dir", "lib"), ("work_dir", "work"), ("output_dir", "output"))}
        self.cfg.update(root=str(self.base), tools_dir=str(Path(TOOLS).resolve()),
                        reference_taxid_map=str(self.base / "mapping.tsv"),
                        masking=True, min_free_gib=0, prepare_target_gib=0.0001,
                        prepare_max_files=1, prepare_threads=2, build_threads=2,
                        build_cores=2, build_memory_gib=16, build_hours=120,
                        build_soft_hours=110, checkpoint_headroom_hours=10,
                        checkpoint_hours=4, build_queue="local-test")
        self.environment = patch.dict(os.environ, OMP_NUM_THREADS="2", K2_DB_READ_THREADS="2")
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.make_references()
        hpc_reference.metadata(self.cfg, source=self.base / "source", download=False)
        plan = prepare.plan(self.cfg)
        self.assertEqual(len(plan["tasks"]), 10)
        for task in plan["tasks"]:
            prepare.prepare_task(self.cfg, task["task_id"])
        self.prepared = prepare.finalize_preparation(self.cfg)
        self.assertEqual(set(self.prepared["groups"]), set(build.CATEGORIES))
        self.assertTrue(all(len(group["inputs"]) == 2 for group in self.prepared["groups"].values()))

    def make_references(self):
        source = self.base / "source"
        source.mkdir()
        nodes = [(1, 1, "no rank"), (2, 1, "superkingdom")]
        mapping = ["path\ttaxid\n"]
        rng = random.Random(837)
        for number, category in enumerate(build.CATEGORIES, 10):
            nodes.append((number, 2, "species"))
            for part in range(2):
                path = self.base / "raw" / category / f"reference{number}_{part}.fna"
                path.parent.mkdir(parents=True, exist_ok=True)
                seq = "".join(rng.choices("ACGT", k=25000))
                path.write_text(f">reference{number}_{part}\n" + "\n".join(seq[i:i+80] for i in range(0, len(seq), 80)) + "\n")
                mapping.append(f"{category}/{path.name}\t{number}\n")
        (self.base / "mapping.tsv").write_text("".join(mapping))
        (source / "nodes.dmp").write_text("".join(
            f"{t}\t|\t{p}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            for t, p, rank in nodes))
        (source / "names.dmp").write_text("".join(
            f"{t}\t|\ttaxon {t}\t|\t\t|\tscientific name\t|\n" for t, _, _ in nodes))
        (source / "merged.dmp").write_text("")

    def checkpoint(self, category="microbes"):
        folder = build.index_root(self.cfg) / category
        pointer = (folder / "checkpoints/CURRENT").read_text().split("\t")[0]
        return json.loads((folder / "checkpoints" / pointer / "receipt.json").read_text())

    def test_pause_resume_budget_change_and_complete_reuse(self):
        paused = {**self.cfg, "native_stop_after_files": 1}
        self.assertEqual(build.build_category(paused, "microbes"), 75)
        folder = build.index_root(self.cfg) / "microbes"
        self.assertFalse((folder / "database.json").exists())
        identity = json.loads((folder / "context.json").read_text())["identity"]
        checkpoint = self.checkpoint()
        self.assertEqual(checkpoint["completed_files"], 1)
        estimate = (folder / "capacity.json").read_bytes()
        changed_resources = {**self.cfg, "build_queue": "different", "build_memory_gib": 20,
                             "build_hours": 112, "build_soft_hours": 100, "checkpoint_hours": 3}
        self.assertEqual(build.build_category(changed_resources, "microbes"), 0)
        self.assertEqual(json.loads((folder / "context.json").read_text())["identity"], identity)
        self.assertEqual((folder / "capacity.json").read_bytes(), estimate)
        self.assertTrue(self.checkpoint()["complete"])
        log = (folder / "build.log").read_text()
        self.assertIn("METAG_RESUME files=1/2", log)
        self.assertEqual(log.count("METAG_FILE_BEGIN 1/2"), 1)
        before = (folder / "build.log").read_bytes()
        self.assertEqual(build.build_category(changed_resources, "microbes"), 0)
        self.assertEqual((folder / "build.log").read_bytes(), before)

    def test_thread_change_and_completed_input_mutation_preserve_checkpoint(self):
        self.assertEqual(build.build_category({**self.cfg, "native_stop_after_files": 1}, "microbes"), 75)
        folder = build.index_root(self.cfg) / "microbes"
        previous = (folder / "checkpoints/CURRENT").read_bytes()
        with self.assertRaisesRegex(RuntimeError, "thread identity changed"):
            build.build_category({**self.cfg, "build_threads": 1}, "microbes")
        self.assertEqual((folder / "checkpoints/CURRENT").read_bytes(), previous)
        path = Path(self.prepared["groups"]["microbes"]["inputs"][0]["path"])
        content = path.read_bytes()
        offset = content.index(b"\n") + 1
        replacement = b"A" if content[offset:offset+1] != b"A" else b"C"
        path.write_bytes(content[:offset] + replacement + content[offset+1:])
        with self.assertRaisesRegex(RuntimeError, "SHA256|changed|checksum"):
            build.build_category(self.cfg, "microbes")
        self.assertEqual((folder / "checkpoints/CURRENT").read_bytes(), previous)

    def test_publish_requires_all_five_and_shared_native_taxonomy(self):
        self.assertEqual(build.build_category(self.cfg, "microbes"), 0)
        with self.assertRaises((FileNotFoundError, RuntimeError)):
            build.publish(self.cfg)
        for category in build.CATEGORIES[1:]:
            self.assertEqual(build.build_category(self.cfg, category), 0)
        published = build.publish(self.cfg)
        self.assertEqual(set(published["context"]["categories"]), set(build.CATEGORIES))
        self.assertEqual(len(published["hashes"]), 5)
        self.assertEqual(published["context"]["absent_categories"], [])
        folder = build.index_root(self.cfg)
        for name in ("taxo.k2d", "opts.k2d"):
            self.assertEqual(len({sha256(folder / category / name) for category in build.CATEGORIES}), 1)
            self.assertEqual(sha256(folder / name), sha256(folder / "microbes" / name))
        self.assertEqual(build.database_snapshot(self.cfg, full=True)["identity"], published["identity"])
        self.assertEqual(build.publish(self.cfg)["identity"], published["identity"])


if __name__ == "__main__":
    unittest.main()
