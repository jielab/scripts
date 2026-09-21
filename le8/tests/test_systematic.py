"""Aggregate audits must retain nulls, missing modules and source tables."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("systematic", Path(__file__).parents[1]/"f/c5_systematic.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class SystematicAudit(unittest.TestCase):
    def test_missing_tests_stay_in_family(self):
        q = audit.bh([.001, .01, np.nan], 100)
        np.testing.assert_allclose(q[:2], [.1, .5])
        self.assertTrue(np.isnan(q[2]))

    def test_no_c4_is_unavailable_not_negative_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = audit.aggregate(Path(tmp))
            self.assertEqual(len(out["inventory"]), 12)
            self.assertTrue((out["inventory"].status == "unavailable").all())
            self.assertTrue(out["all_c4_contrasts"].empty)

    def test_full_family_and_null_pgs_retained_without_source_mutation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); p = root/"prot/c1_correlate"; p.mkdir(parents=True)
            d = pd.DataFrame(dict(term=["A","B","C","D"], beta=[1.,.2,0.,0.],
                **{"std.error":[.1]*4,"p.value":[.001,.1,.8,.9]},
                FDR=[.004,.2,.9,.9],N_total=[1000]*4,N_event=[40]*4))
            d.to_csv(p/"pwas_incident_adj2.csv", index=False)
            g=d.copy();g.FDR=1;g.to_csv(p/"pwas_pgs_incident_full_genetic.csv", index=False)
            lm=d.iloc[:1].copy();lm["landmark_years"]=5;lm.to_csv(p/"pwas_incident_landmark_adj2.csv", index=False)
            before={f.name:f.read_bytes() for f in p.iterdir()}
            out=audit.aggregate(root)
            self.assertEqual(out["c1_summary"].query("analysis=='biomarker_PGS'").FDR05.iloc[0],0)
            self.assertEqual(len(out["measured_PGS"]),4)
            self.assertAlmostEqual(out["landmark_family"].FDR_assay_landmark_family.iloc[0],.004)
            self.assertTrue(out["landmark_family"].subset_selected_before_landmark.iloc[0])
            self.assertEqual(before,{f.name:f.read_bytes() for f in p.iterdir()})


if __name__ == "__main__":
    unittest.main()
