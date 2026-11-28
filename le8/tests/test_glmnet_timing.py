"""Regression checks for the single-path solver and independent timing publication."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

F = Path(__file__).resolve().parents[1] / "f"
sys.path.insert(0, str(F))
import dx_compare
from outputs import default_cohort_path


class TimingPublication(unittest.TestCase):
    def test_failed_parent_keeps_independent_studies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            aod = root/"cad/aod"
            (aod/"main/results").mkdir(parents=True)
            pd.DataFrame([dict(dx="any1",endpoint_sha256="matched")]).to_csv(aod/"main/results/dx_counts.tsv",sep="\t",index=False)
            base = root/"scope"
            for sub in ["c1_correlate","logs","qc"]:
                (base/sub).mkdir(parents=True)
            (base/"qc/dx.json").write_text(json.dumps(dict(trait="cad",dx="any1",endpoint_sha256="matched")))
            (base/"logs/analysis.status.json").write_text('{"status":"FAILED"}')
            labels = dict(trait="cad",layer="prot",dx="any1",endpoint_sha256="matched")
            with pd.ExcelWriter(base/"c1_correlate/timing.xlsx") as book:
                pd.DataFrame([{**labels,"analysis":a} for a in ("controlled_training","saved_C1_common_test")]).to_excel(book,sheet_name="timing_paired",index=False)
                pd.DataFrame([{**labels,"model":"omics"}]).to_excel(book,sheet_name="timing_scores",index=False)
                pd.DataFrame([{**labels,"endpoint_sha256":"changed"}]).to_excel(book,sheet_name="timing_medication",index=False)
            with patch.object(dx_compare,"scope_for_read",return_value=base):
                out = dx_compare.report(root,["cad"],["prot"],["any1"])
            self.assertEqual(list(out["timing_paired"].analysis),["controlled_training"])
            self.assertEqual(len(out["timing_scores"]),1)
            self.assertTrue(out["timing_medication"].empty)
            self.assertTrue(out["performance"].empty)
            roster = default_cohort_path(root)
            self.assertEqual(roster.parent,root/"tmp/freeze")
            old = root/"logs/freeze"/roster.name
            old.parent.mkdir(parents=True)
            old.touch()
            self.assertEqual(default_cohort_path(root),roster)

    def test_real_glmnet_paths(self):
        r = os.environ.get("R_BIN") or shutil.which("Rscript")
        if not r or subprocess.run([r,"-e",'quit(status=as.integer(!all(sapply(c("glmnet","jsonlite"),requireNamespace,quietly=TRUE))))'],
                                   capture_output=True).returncode:
            self.skipTest("R glmnet and jsonlite are required")
        rng = np.random.default_rng(65)
        x = rng.normal(size=(620,50))
        x[:,20:] = x[:,:30] + .01*rng.normal(size=(620,30))
        signal = 1.5*x[:,0]+.5*x[:,1]
        truth = (rng.random(len(x))<1/(1+np.exp(-signal))).astype(int)
        roles = np.repeat(["build","tune_model","test"],[400,150,70])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            x.astype("<f8").tofile(root/"X.bin")
            for label,budget in [("full",5000),("prefix",150),("boundary",12)]:
                out = root/label
                out.mkdir()
                y = truth.copy()
                if label=="prefix":
                    y[400:550] = x[400:550,0]<0
                if label=="boundary":
                    y[400:550] = signal[400:550]>0
                y[550:] = 0
                pd.DataFrame(dict(role=roles,y=y,weight=(roles!="test").astype(int))).to_csv(root/"roles.tsv",sep="\t",index=False)
                (root/"config.json").write_text(json.dumps(dict(seed=3,n=len(x),p=x.shape[1],maxit=budget,feature_names=list(map(str,range(50))))))
                run = subprocess.run([r,str(F/"c1_unified_glmnet.R"),str(root),str(out)],text=True,capture_output=True)
                self.assertEqual(run.returncode,0,run.stderr)
                audit = json.loads((out/"glmnet_fit.json").read_text())
                self.assertEqual(audit["fit_attempts"],1)
                self.assertLessEqual(audit["coordinate_passes"],budget+1)
                self.assertGreaterEqual(audit["solver_seconds"],0)
                self.assertEqual(audit["path_status"]=="solver_completed",label=="full")
                self.assertEqual((out/"raw_probability_LOCAL_ONLY.bin").exists(),label!="boundary")
                if label=="prefix":
                    self.assertGreaterEqual(audit["lambda_count"]-audit["selected_index"],5)
                if label=="boundary":
                    self.assertEqual(audit["status"],"unavailable_lambda_selection")


if __name__ == "__main__":
    unittest.main()
