import copy
import json
from pathlib import Path
import tempfile
import unittest
import joblib
import pandas as pd
from abm_core.model import file_sha256, write_json
from migrate_abm_pooled import (project_identity, write_predictor_cache, verify_files,
                                COMPARATOR, OLD_COMPARATOR)


class CacheMigrationTests(unittest.TestCase):
    def test_dependency_changes_are_rejected_before_cache_projection(self):
        old = dict(ancestry="old-labels", sources={COMPARATOR:OLD_COMPARATOR,
            "grid/f/3.disco.R":"upstream", "grid/f/abm_core/predictors.py":"math"},
            supervised_labels="labels", roles="roles", side_values="features")
        sources = {COMPARATOR:"new-pooled", "grid/f/abm_core/predictors.py":"math"}
        result = project_identity(old, sources)
        self.assertNotIn("ancestry", result)
        self.assertEqual(result["supervised_labels"], old["supervised_labels"])
        self.assertEqual(result["schema_version"], "abm.predictors.pooled.1")
        with self.assertRaisesRegex(ValueError, "dependency changed"):
            project_identity(old, dict(sources, **{"grid/f/abm_core/predictors.py":"changed"}))
        other = copy.deepcopy(old)
        other["sources"][COMPARATOR] = "unknown"
        with self.assertRaisesRegex(ValueError, "audited"):
            project_identity(other, sources)

    def test_cache_migration_keeps_pooled_predictions_and_checks_written_hashes(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            old = root / "old/_LOCAL_ONLY"
            (old / "_LOCAL_ONLY").mkdir(parents=True)
            bundle = dict(base="fixed-baseline", candidate="fixed-candidate",
                          comparators={"Clinical":"retained", "retired_by_ancestry":"unused"})
            joblib.dump(bundle, old / "_LOCAL_ONLY/predictor_bundle.joblib")
            predictions = pd.DataFrame(dict(id=["01","02"], group_id=["a","b"], role="test",
                p_base=[.1,.2], p_candidate=[.3,.4], Clinical=[.5,.6], retired_by_ancestry=[.7,.8]))
            predictions.to_parquet(old / "_LOCAL_ONLY/predictor_predictions.parquet", index=False)
            pd.DataFrame(dict(method=["Clinical","retired_by_ancestry"], status="completed")).to_csv(
                old / "comparator_status.tsv", sep="\t", index=False)
            write_json(old / "predictor_identity.json", {"legacy":True})
            names = ["_LOCAL_ONLY/predictor_bundle.joblib", "_LOCAL_ONLY/predictor_predictions.parquet",
                     "comparator_status.tsv", "predictor_identity.json"]
            saved = dict(executed_device="cpu", files={n:file_sha256(old/n) for n in names})
            write_json(old / "COMPLETE.json", saved)
            new = root / "new/_LOCAL_ONLY"
            identity = {"schema_version":"abm.predictors.pooled.1"}
            write_predictor_cache(old, new, saved, identity)
            migrated = json.loads((new / "COMPLETE.json").read_text())
            verify_files(new, migrated["files"])
            verify_files(old, saved["files"])
            observed = pd.read_parquet(new / "_LOCAL_ONLY/predictor_predictions.parquet")
            pd.testing.assert_frame_equal(observed, predictions.drop(columns="retired_by_ancestry"))
            observed_bundle = joblib.load(new / "_LOCAL_ONLY/predictor_bundle.joblib")
            self.assertEqual(observed_bundle["base"], bundle["base"])
            self.assertEqual(observed_bundle["candidate"], bundle["candidate"])
            self.assertEqual(set(observed_bundle["comparators"]), {"Clinical"})
            # Repeating the migration verifies the destination, never trusts a marker alone.
            write_predictor_cache(old, new, saved, identity)
            (new / "comparator_status.tsv").write_text("corrupted")
            with self.assertRaisesRegex(ValueError, "artifact changed"):
                write_predictor_cache(old, new, saved, identity)
