"""Real reference streams: adoption, whole-record separation and restart guards."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_reference as old
import partitioned_reference as new
from hpc_common import atomic_json, sha256


SOFTWARE = {"source_commit": old.SOURCE_COMMIT, "source_tag": "2.17.2",
            "files": {name: "a" * 64 for name in ("k2mask", "estimate_capacity", "build_db", "classify")}}


class PartitionedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cfg = {key: str(self.root / name) for key, name in (
            ("reference_raw", "raw"), ("metadata_dir", "metadata"),
            ("library_dir", "lib"), ("tools_dir", "tools"))}
        self.cfg.update(root=str(self.root), masking=False, min_free_gib=0,
                        prepare_target_gib=1, prepare_max_files=1000, prepare_threads=2,
                        reference_taxid_map=str(self.root / "mapping.tsv"))
        self.software = patch.object(old, "_software", return_value=SOFTWARE)
        self.software.start()
        self.addCleanup(self.software.stop)

    def fixture(self, categories=("microbes", "viruses", "plants", "vertebrates", "invertebrates"), *, old_cache=True, task_files=1000):
        source = self.root / "source"
        source.mkdir()
        nodes = [(1, 1, "no rank"), (2, 1, "superkingdom")]
        for index, category in enumerate(categories, 10):
            path = self.root / "raw" / category / f"reference{index}.fna"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f">old id {index}/a\nACGTACGT\nACGT\n>old id {index}/b\nNNGGCCAATT\n")
            nodes.append((index, 2, "species"))
        (source / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{rank}\t|\n" for t, p, rank in nodes))
        (source / "names.dmp").write_text("".join(f"{t}\t|\ttaxon {t}\t|\t\t|\tscientific name\t|\n" for t, _, _ in nodes))
        (source / "merged.dmp").write_text("")
        (self.root / "mapping.tsv").write_text("path\ttaxid\n" + "".join(
            f"{category}/reference{index}.fna\t{index}\n" for index, category in enumerate(categories, 10)))
        old.metadata(self.cfg, source=source, download=False)
        if old_cache:
            old_cfg = {**self.cfg, "prepare_max_files": task_files}
            planned = old.plan_references(old_cfg)
            for task in planned["tasks"]:
                old.prepare_task(old_cfg, task["task_id"])
            return planned

    def prepared(self, cfg=None):
        cfg = cfg or self.cfg
        planned = new.plan(cfg)
        for task in planned["tasks"]:
            new.prepare_task(cfg, task["task_id"])
        return new.finalize_preparation(cfg)

    def test_mixed_old_task_splits_whole_records_without_masking(self):
        old_plan = self.fixture()
        self.assertEqual(len(old_plan["tasks"]), 1)
        with patch.object(old, "_prepare_native", side_effect=AssertionError("must reuse mask output")):
            catalog = self.prepared()
        self.assertEqual(set(catalog["groups"]), set(new.CATEGORIES))
        global_map = Path(catalog["seqid_map"]).read_text()
        self.assertEqual(len(global_map.splitlines()), 10)
        seen = set()
        for category, group in catalog["groups"].items():
            self.assertEqual((group["records"], group["bases"], group["longest_record"]), (2, 22, 12))
            self.assertEqual(len(group["inputs"]), 1)
            row = group["inputs"][0]
            path = Path(row["path"])
            self.assertEqual(path.name, category + ".fna")
            self.assertEqual(sha256(path), row["sha256"])
            lines = path.read_text().splitlines()
            headers = [line[1:] for line in lines if line.startswith(">")]
            self.assertEqual(len(headers), 2)
            self.assertFalse(seen & set(headers))
            seen.update(headers)
            self.assertIn("ACGTACGT\nACGT\n", path.read_text())
            self.assertEqual(Path(group["inputs_tsv"]).read_text(), "\t".join(str(row[k]) for k in ("path", "size", "sha256", "records", "bases")) + "\n")
        self.assertEqual(seen, {line.split("\t")[0] for line in global_map.splitlines()})
        self.assertEqual(new.catalog_snapshot(self.cfg)["identity"], catalog["identity"])

    def test_whole_category_old_tasks_are_referenced_and_not_copied(self):
        self.fixture(task_files=1)
        catalog = self.prepared()
        for group in catalog["groups"].values():
            self.assertIn("/lib/prepared/task", group["inputs"][0]["path"])
            self.assertNotIn("/parallel/", group["inputs"][0]["path"])
        for folder in (self.root / "lib/parallel/prepared").glob("task*"):
            if folder.is_dir():
                receipt = json.loads((folder / "receipt.json").read_text())
                self.assertEqual(receipt["execution"]["historical_integrity"], new.BASELINE_ONLY)
                self.assertFalse((folder / "sequences.fna").exists())

    def test_plan_adopts_copied_legacy_and_ignores_old_absolute_paths(self):
        self.fixture(task_files=1)
        moved = self.root / "moved"
        shutil.copytree(self.root / "raw", moved / "raw")
        shutil.copytree(self.root / "metadata", moved / "meta")
        shutil.copytree(self.root / "lib", moved / "oldlib")
        shutil.rmtree(self.root / "raw")
        shutil.rmtree(self.root / "metadata")
        shutil.rmtree(self.root / "lib")
        (self.root / "mapping.tsv").unlink()
        cfg = {**self.cfg, "reference_raw": str(moved / "raw"), "metadata_dir": str(moved / "freshmeta"),
               "library_dir": str(moved / "newlib"), "legacy_library_dir": str(moved / "oldlib"),
               "legacy_metadata_dir": str(moved / "meta")}
        with patch.object(old, "_prepare_native", side_effect=AssertionError("adoption must not remask")):
            catalog = self.prepared(cfg)
        self.assertTrue(all(str(moved) in row["path"] for group in catalog["groups"].values() for row in group["inputs"]))
        self.assertEqual(sum(group["records"] for group in catalog["groups"].values()), 10)

    def test_complete_parallel_cache_is_portable_with_stable_content_identity(self):
        self.fixture(task_files=1)
        initial = self.prepared()
        moved = self.root / "another"
        for name in ("raw", "metadata", "lib"):
            shutil.copytree(self.root / name, moved / name)
        cfg = {**self.cfg, "reference_raw": str(moved / "raw"), "metadata_dir": str(moved / "metadata"), "library_dir": str(moved / "lib")}
        for name in ("raw", "metadata", "lib"):
            shutil.rmtree(self.root / name)
        self.assertEqual(new.plan(cfg)["identity"], new.load_plan(cfg)["identity"])
        with self.assertRaisesRegex(RuntimeError, "relocated"):
            new.catalog_snapshot(cfg)
        final = new.finalize_preparation(cfg)
        self.assertEqual(final["identity"], initial["identity"])
        self.assertEqual(new.catalog_snapshot(cfg)["identity"], final["identity"])

    def test_completed_receipt_reuses_stat_guard_without_rehashing_prepared(self):
        self.fixture(task_files=1)
        self.prepared()
        actual = new.sha256
        calls = []
        def tracked(path):
            calls.append(str(path))
            if str(path).endswith(("sequences.fna", "seqid2taxid.map")):
                raise AssertionError("completed prepared data should not be rehashed")
            return actual(path)
        with patch.object(new, "sha256", side_effect=tracked):
            receipt = new.prepare_task(self.cfg, 1)
        self.assertEqual(receipt["status"], "complete")
        self.assertTrue(calls)

    def test_same_size_prepared_mutation_after_seal_is_rejected(self):
        self.fixture(task_files=1)
        self.prepared()
        path = self.root / "lib/prepared/task000001/sequences.fna"
        data = path.read_bytes()
        path.write_bytes(data.replace(b"ACGTACGT", b"TCGTACGT", 1))
        with self.assertRaisesRegex(RuntimeError, "SHA256"):
            new.prepare_task(self.cfg, 1)

    def test_bad_legacy_record_counter_is_rejected_at_adoption(self):
        self.fixture()
        new.plan(self.cfg)
        path = self.root / "lib/prepared/task000001/sequences.fna"
        path.write_bytes(path.read_bytes().replace(b"_r0000000002", b"_r0000000003", 1))
        with self.assertRaisesRegex(RuntimeError, "record counter"):
            new.prepare_task(self.cfg, 1)
        self.assertFalse((self.root / "lib/parallel/prepared/task000001/receipt.json").exists())

    def test_partial_owned_stage_restarts_only_one_task(self):
        self.fixture(task_files=1)
        planned = new.plan(self.cfg)
        first = new.prepare_task(self.cfg, 1)
        base = self.root / "lib/parallel/prepared"
        partial = base / ".task000002.partial"
        partial.mkdir()
        atomic_json(partial / "owner.json", {"schema": new.SCHEMA, "kind": "parallel-reference-task",
                    "plan_identity": planned["identity"], "task_id": 2})
        (partial / "junk").write_text("interrupted")
        new.prepare_task(self.cfg, 2)
        self.assertFalse(partial.exists())
        self.assertEqual(new.prepare_task(self.cfg, 1)["identity"], first["identity"])

    def test_foreign_partial_is_not_removed(self):
        self.fixture()
        new.plan(self.cfg)
        partial = self.root / "lib/parallel/prepared/.task000001.partial"
        partial.mkdir(parents=True)
        (partial / "foreign.txt").write_text("keep")
        with self.assertRaisesRegex(RuntimeError, "unowned"):
            new.prepare_task(self.cfg, 1)
        self.assertEqual((partial / "foreign.txt").read_text(), "keep")

    def test_fresh_reference_plan_preserves_global_ids_across_categories(self):
        self.fixture(old_cache=False)
        catalog = self.prepared()
        self.assertEqual(len(new.load_plan(self.cfg)["tasks"]), 5)
        entries = Path(catalog["seqid_map"]).read_text().splitlines()
        ids = [int(line.split("_f")[1].split("_")[0]) for line in entries]
        self.assertEqual(ids, [0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        self.assertEqual(new.catalog_snapshot(self.cfg, full=True)["identity"], catalog["identity"])

    def test_missing_legacy_receipt_prepares_only_missing_task(self):
        self.fixture(task_files=1)
        shutil.rmtree(self.root / "lib/prepared/task000003")
        planned = new.plan(self.cfg)
        self.assertEqual([t["mode"] for t in planned["tasks"]], ["adopt", "adopt", "prepare", "adopt", "adopt"])
        with patch.object(old, "_prepare_native", wraps=old._prepare_native) as preparation:
            for task in planned["tasks"]:
                new.prepare_task(self.cfg, task["task_id"])
        self.assertEqual(preparation.call_count, 1)
        self.assertEqual(sum(g["records"] for g in new.finalize_preparation(self.cfg)["groups"].values()), 10)

    def test_unknown_category_requires_explicit_mapping(self):
        self.fixture(categories=("other",), old_cache=False)
        with self.assertRaisesRegex(RuntimeError, "Unknown reference category"):
            new.plan(self.cfg)
        cfg = {**self.cfg, "reference_category_map": {"other": "microbes"}}
        self.assertEqual(new.plan(cfg)["categories"], ["microbes"])

    def test_singular_category_aliases_are_supported(self):
        self.fixture(categories=("plant", "viral", "vertebrate", "invertebrate", "microbe"), old_cache=False)
        self.assertEqual(set(new.plan(self.cfg)["categories"]), set(new.CATEGORIES))

    def test_modified_manifest_or_taxid_cannot_be_adopted(self):
        self.fixture()
        path = self.root / "lib/reference_manifests/task000001.json"
        content = path.read_text()
        path.write_text(content.replace('"taxid": 14', '"taxid": 15'))
        with self.assertRaisesRegex(RuntimeError, "SHA256"):
            new.plan(self.cfg)

    def test_global_taxonomy_or_native_list_mutation_is_rejected(self):
        self.fixture()
        catalog = self.prepared()
        path = Path(catalog["groups"]["plants"]["inputs_tsv"])
        original = path.read_text()
        path.write_text(original.replace("\t2\t22", "\t3\t22"))
        with self.assertRaisesRegex(RuntimeError, "input list changed"):
            new.catalog_snapshot(self.cfg)
        path.write_text(original)
        taxonomy = Path(catalog["taxonomy_dir"]) / "names.dmp"
        taxonomy.write_bytes(taxonomy.read_bytes().replace(b"taxon", b"TAXON", 1))
        with self.assertRaisesRegex(RuntimeError, "SHA256"):
            new.catalog_snapshot(self.cfg)

    def test_raw_inventory_change_blocks_reuse(self):
        self.fixture()
        new.plan(self.cfg)
        path = self.root / "raw/microbes/reference10.fna"
        path.write_text(path.read_text() + "ACTG\n")
        with self.assertRaisesRegex(RuntimeError, "inventory"):
            new.plan(self.cfg)

    def test_incoming_filepart_blocks_inventory_even_with_finished_copy(self):
        raw = self.root / "raw/microbes"
        raw.mkdir(parents=True)
        (raw / "reference.fna.gz").write_bytes(b"existing")
        (raw / "reference.fna.gz.filepart").write_bytes(b"uploading")
        with self.assertRaisesRegex(RuntimeError, "Unfinished reference download/file.*filepart"):
            old._inventory(self.cfg)

    def test_rsync_staging_fasta_is_not_accepted_as_another_reference(self):
        partial = self.root / "raw/microbes/.rsync-partial-qiming-20261104"
        partial.mkdir(parents=True)
        (partial / "GCF_000001.1_genomic.fna.gz").write_bytes(b"truncated")
        with self.assertRaisesRegex(RuntimeError, "Rsync staging directory"):
            old._inventory(self.cfg)


if __name__ == "__main__":
    unittest.main()
