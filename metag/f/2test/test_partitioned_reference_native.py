"""Optional real CPU k2mask integration; no mocked native receipts.

METAG_NATIVE_FIXTURE must name a tiny fixture created with build_native_fixture.py.
Only that fixture's original installed CPU tools and two reference inputs are
read. All preparation/adoption outputs go to a new temporary directory.
"""
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_reference
import partitioned_reference
from hpc_common import sha256


@unittest.skipUnless(os.environ.get("METAG_NATIVE_FIXTURE"), "set METAG_NATIVE_FIXTURE for the real installed CPU tools")
class NativePartitionedPreparationTests(unittest.TestCase):
    def test_real_masked_fresh_preparation_and_old_cache_adoption(self):
        fixture = Path(os.environ["METAG_NATIVE_FIXTURE"]).resolve()
        fixture_cfg = json.loads((fixture / "config.json").read_text())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw"
            for name, taxid in (("microbes", 10), ("viruses", 20)):
                target = raw / name / f"species{taxid}.fna"
                target.parent.mkdir(parents=True)
                shutil.copyfile(fixture / "raw" / f"species{taxid}.fna", target)
            mapping = root / "map.tsv"
            mapping.write_text("path\ttaxid\nmicrobes/species10.fna\t10\nviruses/species20.fna\t20\n")
            cfg = {**fixture_cfg, "root": str(root), "reference_raw": str(raw),
                   "metadata_dir": str(root / "metadata"), "library_dir": str(root / "newlib"),
                   "tools_dir": fixture_cfg["tools_dir"], "reference_taxid_map": str(mapping),
                   "masking": True, "min_free_gib": 0, "prepare_threads": 2,
                   "prepare_target_gib": 1, "prepare_max_files": 1000}
            hpc_reference.metadata(cfg, source=fixture / "source", download=False)
            fresh_plan = partitioned_reference.plan(cfg)
            self.assertEqual(fresh_plan["categories"], ["microbes", "viruses"])
            for task in fresh_plan["tasks"]:
                receipt = partitioned_reference.prepare_task(cfg, task["task_id"])
                self.assertEqual(receipt["execution"]["mode"], "prepare")
            fresh = partitioned_reference.finalize_preparation(cfg)
            self.assertEqual(sum(group["bases"] for group in fresh["groups"].values()), 200000)
            self.assertEqual(sum(group["records"] for group in fresh["groups"].values()), 2)
            old_cfg = {**cfg, "library_dir": str(root / "oldlib")}
            old_plan = hpc_reference.plan_references(old_cfg)
            self.assertEqual(len(old_plan["tasks"]), 1)
            hpc_reference.prepare_task(old_cfg, 1)
            adopt_cfg = {**cfg, "library_dir": str(root / "adoptedlib"),
                         "legacy_library_dir": str(root / "oldlib"), "legacy_metadata_dir": str(root / "metadata")}
            adopted_plan = partitioned_reference.plan(adopt_cfg)
            self.assertEqual([task["mode"] for task in adopted_plan["tasks"]], ["adopt"])
            adopted_receipt = partitioned_reference.prepare_task(adopt_cfg, 1)
            self.assertFalse(adopted_receipt["execution"]["masking_repeated"])
            adopted = partitioned_reference.finalize_preparation(adopt_cfg)
            self.assertEqual(Path(fresh["seqid_map"]).read_bytes(), Path(adopted["seqid_map"]).read_bytes())
            for category in ("microbes", "viruses"):
                left = fresh["groups"][category]["inputs"][0]
                right = adopted["groups"][category]["inputs"][0]
                self.assertEqual(left["sha256"], right["sha256"])
                self.assertEqual(sha256(right["path"]), right["sha256"])
            self.assertEqual(partitioned_reference.catalog_snapshot(adopt_cfg, full=True)["identity"], adopted["identity"])


if __name__ == "__main__":
    unittest.main()
