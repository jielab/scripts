import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from abm_core.model import file_sha256
from kinship import active, correct_benchmark, correct_grid
from prepare_kinship import activate_correction
from prsformer import read_inputs, split_registry, validate_family_splits, validate_prediction_holdout


def fixture(root):
    cache = root / "prsformer/cache"
    cache.mkdir(parents=True)
    data = pd.DataFrame(dict(eid=[str(i) for i in range(104)], group=[str(i) for i in range(104)],
                             target="ALL", split=["train"] * 80 + ["validation"] * 20 + ["test"] * 4,
                             height=np.arange(104, dtype=float)))
    data.to_csv(cache / "data.tsv.gz", sep="\t", index=False)
    np.save(cache / "genotypes.npy", np.zeros((104, 4), np.float16))
    pd.DataFrame(dict(CHR=1, BP=np.arange(1, 5), SNP=[f"rs{x}" for x in range(4)],
                      REF="A", ALT="C")).to_csv(cache / "variants.tsv.gz", sep="\t", index=False)
    out = root / "local/kinship"
    out.mkdir(parents=True)
    groups = data.group.copy()
    groups.iloc[0] = groups.iloc[100]  # A development relative of a test person.
    groups.iloc[3] = groups.iloc[2]  # A family retained wholly in development.
    pd.DataFrame(dict(eid=data.eid, group_id=groups)).to_csv(out / "groups.tsv", sep="\t", index=False)
    pd.DataFrame(dict(eid=data.eid, group_id=groups, is_test=data.split.eq("test"),
                     eligible=[False] + [True] * 103)).to_csv(out / "eligible.tsv", sep="\t", index=False)
    source = out / "source.dat"
    source.write_text("ID1 ID2 Kinship\n0 100 0.1\n2 3 0.1\n")
    evidence = dict(status="verified", source=str(out / "groups.tsv"), sha256=file_sha256(out / "groups.tsv"))
    (out / "group_verification.json").write_text(json.dumps(evidence))
    activate_correction(root, out, source, out / "groups.tsv", out / "eligible.tsv")
    return cache, data, out


class KinshipOverlayTests(unittest.TestCase):
    def test_grid_excludes_relatives_and_keeps_exact_test_roster(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, data, out = fixture(root)
            frame = data.rename(columns={"eid": "id", "group": "group_id"})
            kept, metadata = correct_grid(frame, {}, root)
            self.assertEqual(set(kept.loc[kept.split.eq("test"), "id"]), set(data.loc[data.split.eq("test"), "eid"]))
            self.assertNotIn("0", set(kept.id))
            self.assertEqual(metadata["excluded_development_relatives"], 1)
            self.assertEqual(metadata["group_verification"]["status"], "verified")
            self.assertEqual(kept.set_index("id").loc["2", "group_id"], kept.set_index("id").loc["3", "group_id"])

    def test_benchmark_keeps_genotype_positions_but_never_trains_excluded_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache, original, _ = fixture(root)
            args = SimpleNamespace(data=str(cache / "data.tsv.gz"), variants=str(cache / "variants.tsv.gz"),
                                   genotypes=str(cache / "genotypes.npy"), seed=20260904)
            data, labels, identity, _ = read_inputs(args, ["height"])
            self.assertEqual(data.eid.tolist(), original.eid.tolist())
            self.assertEqual(data.loc[0, "split"], "unused")
            self.assertTrue(data.loc[100:, "split"].eq("test").all())
            self.assertEqual(data.loc[2, "split"], data.loc[3, "split"])
            validate_family_splits(data)
            self.assertEqual(identity["kinship_correction_sha256"], active(root)["signature"])
            self.assertEqual(len(split_registry(data)["development_ids"]), 99)
            with self.assertRaisesRegex(ValueError, "predates"):
                validate_prediction_holdout(data, identity, {"training_identity": {}})

    def test_separate_calibration_can_never_reenable_excluded_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache, original, _ = fixture(root)
            data, _ = correct_benchmark(original, root, 20260904, separate_calibration=True)
            roles = pd.DataFrame(dict(id=data.eid, group_id=data.group,
                                      role=["excluded_development_relative"] + ["build"] * 79 +
                                           ["tune_model"] * 10 + ["prediction_calibration"] * 10 + ["test"] * 4))
            path = root / "roles.tsv"
            roles.to_csv(path, sep="\t", index=False)
            args = SimpleNamespace(data=str(cache / "data.tsv.gz"), variants=str(cache / "variants.tsv.gz"),
                                   genotypes=str(cache / "genotypes.npy"), calibration_roles=str(path))
            d, *_ = read_inputs(args, ["height"])
            self.assertEqual(d.loc[0, "split"], "unused")
            self.assertEqual(d.split.eq("prediction_calibration").sum(), 10)
            roles.loc[0, "role"] = "build"
            roles.to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "Related-family|excluded"):
                read_inputs(args, ["height"])

    def test_changed_evidence_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = fixture(root)
            (out / "source.dat").write_text("changed")
            with self.assertRaisesRegex(ValueError, "kinship input changed"):
                active(root)


if __name__ == "__main__":
    unittest.main()
