"""Migration checks with real files; no UKB data or inference is used."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest

import pandas as pd

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("grid_migration_data", HERE / "grid.data.py")
data = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = data
SPEC.loader.exec_module(data)
common = data.load_module(HERE / "0.common.py", "grid_migration_common")


class HistoricalWeights(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="grid-migration-")
        self.addCleanup(temporary.cleanup)
        self.root = root = Path(temporary.name)
        self.snpinfo = root / "snpinfo"
        self.snpinfo.write_text("CHR\tSNP\tBP\n1\trs1\t101\n")
        self.sources = {pop: root / f"height.{pop}.gz" for pop in data.POPS}
        for source in self.sources.values():
            pd.DataFrame({"SNP": ["rs1"], "A1": ["A"], "A2": ["C"],
                          "BETA": [.1], "SE": [.01], "N": [1000]}).to_csv(source, sep="\t", index=False)
        self.signature = json.dumps({"settings": {"phi": "1e-2", "iterations": 1000,
            "burnin": 500, "thin": 5, "seed": 1, "n_gwas": "", "chromosomes": [1]},
            "files": [common.stamp(p) for p in [*self.sources.values(), self.snpinfo]]})
        records = {}
        self.metadata_paths = {}
        for pop, source in self.sources.items():
            weight = root / f"height.{pop}.chr1.csx.gz"
            pd.DataFrame({"CHR": [1], "BP": [101], "SNP": ["rs1"], "A1": ["A"],
                          "A2": ["C"], "BETA": [.1]}).to_csv(weight, sep="\t", index=False)
            # The original preparation key hashes a previous code version.
            old_key = hashlib.sha256(json.dumps(["height", pop, common.stamp(source),
                common.stamp(self.snpinfo), "historical-preprocessing-code-digest"]).encode()).hexdigest()
            metadata = {"preparation_signature": old_key, "trait": "height", "pop": pop,
                        "input": str(source.resolve()), "n_gwas_used": 1000, "inference_phi": "1e-2"}
            self.metadata_paths[pop] = Path(str(weight) + ".metadata.json")
            self.metadata_paths[pop].write_text(json.dumps(metadata))
            Path(str(weight) + ".signature").write_text(self.signature)
            records.update(common.score_source_records({f"csx.{pop}": weight}, self.signature, [1]))
        (root / "scores/height").mkdir(parents=True)
        score = root / "scores/height/1csx.scores.rds"
        score.write_bytes(b"synthetic score identity; no outcomes")
        self.provenance = score.with_suffix(".provenance.json")
        self.provenance.write_text(json.dumps({"schema": "grid_csx_scores_v1", "models": records,
                                               "score_file_sha256": common.sha256(score)}))
        self.args = SimpleNamespace(gwas_dir=root, snpinfo=self.snpinfo, score_dir=root / "scores", chrs="1")

    def test_rescored_historical_weights_load_without_relabeling_metadata(self):
        before = {p: p.read_bytes() for p in self.metadata_paths.values()}
        weights, sources = data.native_weights(self.args, "height")
        self.assertEqual(len(weights), 1)
        self.assertEqual(weights.filter(like="beta_").iloc[0].tolist(), [.1] * 4)
        self.assertEqual(sum(s.get("preparation_validation") ==
            "historical_preparation_inference_inputs_verified" for s in sources), 4)
        self.assertTrue(all(p.read_bytes() == value for p, value in before.items()))

    def test_historical_weights_still_require_matching_rescored_provenance(self):
        original = self.provenance.read_bytes()
        self.provenance.unlink()
        with self.assertRaisesRegex(ValueError, "Missing CSx score/weight provenance"):
            data.native_weights(self.args, "height")
        records = json.loads(original)
        records["models"]["csx.EUR"]["weights_sha256"] = "different-weights"
        self.provenance.write_text(json.dumps(records))
        with self.assertRaisesRegex(ValueError, "different CSx runs"):
            data.native_weights(self.args, "height")

    def test_changed_gwas_or_reference_is_not_a_code_only_migration(self):
        for source in (self.sources["EUR"], self.snpinfo):
            with self.subTest(source=source.name):
                original = source.read_bytes()
                stat = source.stat()
                try:
                    source.write_bytes(original + b"\n")
                    with self.assertRaisesRegex(ValueError, "preparation mismatch"):
                        data.native_weights(self.args, "height")
                finally:
                    source.write_bytes(original)
                    import os
                    os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))

    def test_unverifiable_historical_metadata_is_rejected(self):
        source = self.sources["EUR"]
        metadata = json.loads(self.metadata_paths["EUR"].read_text())
        for updates in ({"trait": "ldl"}, {"pop": "EAS"}, {"input": "/different/gwas"},
                        {"preparation_signature": ""}):
            with self.subTest(updates=updates), self.assertRaisesRegex(ValueError, "preparation mismatch"):
                data.verify_weight_preparation(common, {**metadata, **updates}, self.signature,
                                               source, self.snpinfo, "height", "EUR")
        for signature in ("unverifiable-old-hash", "{}", "null", "[]"):
            with self.subTest(signature=signature), self.assertRaisesRegex(ValueError, "preparation mismatch"):
                data.verify_weight_preparation(common, metadata, signature, source, self.snpinfo, "height", "EUR")


if __name__ == "__main__":
    unittest.main(verbosity=2)
