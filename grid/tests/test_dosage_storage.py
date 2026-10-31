"""Real dosage -> block -> compact probe -> delete -> consumer round trips."""
import json
import shutil
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import block_inputs as b
import dosage_storage as storage
from abm_core.grid_inputs import locate_blocks
from test_block_producer import fixture
from test_block_dosage import write_pgen


class DosageStorageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.path, self.cfg = fixture(self.root, "ALL.37.bed")
        self.outer = dict(block_config=str(self.path), execution={"device": "cpu"})

    def prepare(self, supplement=False):
        if supplement:
            geno = np.load(self.root / "genotype/genotypes.npy")
            cohort = pd.read_csv(self.root / "genotype/data.tsv.gz", sep="\t", dtype=str)
            write_pgen(self.root / "typ", cohort.eid.tolist(), [("rs2", 150, "C", "T")], [geno[:, 1]])
            np.save(self.root / "genotype/genotypes.npy", geno[:, :1])
            variants = pd.read_csv(self.root / "genotype/variants.tsv.gz", sep="\t")
            variants.iloc[:1].to_csv(self.root / "genotype/variants.tsv.gz", sep="\t", index=False)
            self.cfg["dosage_supplement_dir"] = str(self.root / "typ")
            self.path.write_text(json.dumps(self.cfg))
        _, self.base, _ = b.prepare_trait(self.path, "height", device="cpu")
        arr, _, _, receipt, _ = locate_blocks(self.outer, "height", allow_build=False)
        arr._mmap.close()
        return receipt

    def clean(self, **kwargs):
        return storage.release_dosages(self.path, traits=["height"], analysis_root=self.root, **kwargs)

    def test_release_and_reuse_exact_probe_with_real_supplement(self):
        first = self.prepare(supplement=True)
        frozen = (self.base / "BLOCKS_COMPLETE.json").read_bytes()
        size = self.clean()
        self.assertGreater(size, 0)
        self.assertTrue((self.root / "genotype/genotypes.npy").exists())
        self.assertEqual(self.clean(apply=True), size)
        self.assertFalse((self.root / "genotype/genotypes.npy").exists())
        self.assertEqual(list((self.root / "cache/dosage_supplements").glob("*/genotypes.npy")), [])
        with patch.object(b, "_prepare_trait", side_effect=AssertionError("unexpected rebuild")):
            arr, _, _, reused, _ = locate_blocks(self.outer, "height", allow_build=False)
            arr._mmap.close()
            b.prepare_trait(self.path, "height", device="cpu")
        self.assertEqual(first["adapter_probe"], reused["adapter_probe"])
        self.assertEqual((self.base / "BLOCKS_COMPLETE.json").read_bytes(), frozen)
        self.assertEqual(self.clean(apply=True), 0)

    def test_corrupted_probe_or_block_prevents_deletion(self):
        self.prepare()
        for filename in (storage.PROBE, "block_scores.npy"):
            path = self.base / filename
            before = path.read_bytes()
            path.write_bytes(before + b"corruption")
            with self.assertRaises(ValueError):
                self.clean(apply=True)
            self.assertTrue((self.root / "genotype/genotypes.npy").exists())
            path.write_bytes(before)

    def test_open_genotype_and_active_lock_are_retained(self):
        self.prepare()
        with (self.root / "genotype/genotypes.npy").open("rb"):
            self.assertEqual(self.clean(apply=True), 0)
        with storage.dosage_lock(self.cfg, exclusive=True):
            with self.assertRaises(BlockingIOError):
                self.clean(apply=True)

    def test_missing_probe_and_unfinished_training_are_retained(self):
        self.prepare()
        run = self.root / "prsformer/run"
        run.mkdir(parents=True)
        (run / "training.pt").write_bytes(b"unfinished checkpoint")
        self.assertEqual(self.clean(apply=True), 0)
        (run / "training.pt").unlink()
        (self.base / storage.MARKER).unlink()
        self.assertEqual(self.clean(apply=True), 0)

    def test_configured_missing_trait_preserves_shared_dosage(self):
        self.prepare()
        self.assertEqual(storage.release_dosages(self.path, traits=["height", "ldl"],
                         analysis_root=self.root, apply=True), 0)

    def test_rebuild_preserves_frozen_metadata_and_rejects_changed_dosage(self):
        from prsformer_cleanup import restore_released
        self.prepare()
        cache = self.root / 'genotype'
        work = self.root / 'rebuilt'
        shutil.copytree(cache, work)
        before = {name: (cache / name).read_bytes() for name in
                  ('data.tsv.gz', 'variants.tsv.gz', 'prepare.json')}
        stamp = (cache / 'genotypes.npy').stat().st_mtime_ns
        self.clean(apply=True)
        original = (work / 'genotypes.npy').read_bytes()
        (work / 'genotypes.npy').write_bytes(original + b'corrupt')
        with self.assertRaisesRegex(ValueError, 'Rebuilt dosage differs'):
            restore_released(cache, work)
        self.assertFalse((cache / 'genotypes.npy').exists())
        (work / 'genotypes.npy').write_bytes(original)
        # Recompression need not reproduce gzip timestamps; values must agree.
        frame = pd.read_csv(work / 'data.tsv.gz', sep='\t', dtype=str)
        frame.to_csv(work / 'data.tsv.gz', sep='\t', index=False)
        self.assertTrue(restore_released(cache, work))
        self.assertEqual((cache / 'genotypes.npy').stat().st_mtime_ns, stamp)
        self.assertEqual((cache / 'genotypes.npy').read_bytes(), original)
        for name, value in before.items():
            self.assertEqual((cache / name).read_bytes(), value)


if __name__ == "__main__":
    unittest.main()
