import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from abm_core.contracts import (write_array_manifest, verify_array_binding,
                               write_feature_manifest, verify_feature_binding)


class ContractTests(unittest.TestCase):
    def test_same_length_row_reordering_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "z.npy"
            np.save(path, [[1., 2.], [3., 4.]])
            meta = pd.DataFrame({"id": ["001", "002"], "group_id": ["g1", "g2"]})
            receipt = write_array_manifest(path, meta, geometry_sha256="a" * 64)
            verify_array_binding(path, meta, receipt)
            with self.assertRaisesRegex(ValueError, "row IDs"):
                verify_array_binding(path, meta.iloc[::-1], receipt)

    def test_changed_array_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "z.npy"
            np.save(path, [[1., 2.]])
            meta = pd.DataFrame({"id": ["001"], "group_id": ["g1"]})
            receipt = write_array_manifest(path, meta, geometry_sha256="a" * 64)
            np.save(path, [[7., 2.]])
            with self.assertRaisesRegex(ValueError, "bytes changed"):
                verify_array_binding(path, meta, receipt)

    def test_feature_producer_and_bytes_must_match(self):
        manifest = {"baseline_id": "base", "candidate_id": "candidate", "bank_sha256": "a" * 64,
                    "geometry_sha256": "b" * 64, "predictor_bundle_sha256": "c" * 64,
                    "role_manifest_sha256": "d" * 64}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "features.tsv"
            path.write_text("id\tn_twins\n001\t90\n")
            receipt = write_feature_manifest(path, manifest=manifest)
            verify_feature_binding(path, receipt, manifest=manifest)
            changed = dict(manifest, bank_sha256="f" * 64)
            with self.assertRaisesRegex(ValueError, "bank_sha256"):
                verify_feature_binding(path, receipt, manifest=changed)
            path.write_text("id\tn_twins\n001\t91\n")
            with self.assertRaisesRegex(ValueError, "features changed"):
                verify_feature_binding(path, receipt, manifest=manifest)


if __name__ == "__main__":
    unittest.main()
