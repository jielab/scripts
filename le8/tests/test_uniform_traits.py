"""Default C1 comparisons use the same endpoints and model options across all traits."""
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

F=Path(__file__).resolve().parents[1]/"f"
sys.path.insert(0,str(F))
from workflow_definition import build_jobs, parse_cli, DEFAULT
import le8_scheduler as scheduler


class UniformStudy(unittest.TestCase):
    def test_default_plan_is_symmetric(self):
        traits=["cvd_cad","stroke","cvd_pad","cvd_afib","mi","cvd_hfail"]
        with tempfile.TemporaryDirectory() as temp, patch.dict(os.environ,{"submit":"hpc","WORKFLOW_CORES":"8","SHARED_PYTHON":sys.executable,"PYTHON_BIN":sys.executable}):
            args,extra=parse_cli([DEFAULT,"--analysis-root",temp,"--Y",",".join(traits),"--biom","prot,met"])
            raw=build_jobs(args,extra)
            jobs={j.key:j for j in raw}
            for trait in traits:
                for assay in ("prot","met"):
                    for dx in ("icd","any1"):
                        job=jobs[f"{dx}.{trait}.{assay}.analysis"]
                        steps=json.loads(job.argv[job.argv.index("--steps-json")+1])
                        self.assertIn(trait+"/aod/",job.env["LE8_DX_TABLE"])
                        self.assertEqual(job.env["LE8_DX"],dx)
                        names=[step["key"].split(".",2)[-1] for step in steps]
                        self.assertEqual("c2_cause" in names,dx=="icd")
                        self.assertEqual("c3_coloc" in names,dx=="icd")
                        self.assertEqual("c4_connect" in names,dx=="icd")
                        self.assertEqual("c5_cellulation" in names,dx=="icd")
                        fit=steps[names.index("c1.primary.models")]["argv"]
                        methods=fit[fit.index("--methods")+1].split(",")
                        self.assertTrue({"clinical","glmnet","abm_clinical"}<=set(methods))
                        if dx=="any1":
                            self.assertEqual(set(methods),{"clinical","glmnet","abm_clinical"})
                            self.assertEqual(fit[fit.index("--association")+1],"skip")
            packed=scheduler.join_endpoint_fits(scheduler.join_aod_focus(scheduler.normalized(raw)))
            self.assertEqual(len(packed),3*len(traits))
            for trait in traits:
                self.assertEqual({j["key"] for j in packed if trait in j["key"].split(".")},
                                 {f"aod-focus.{trait}",f"icd-any1.{trait}.prot",f"icd-any1.{trait}.met"})

    def test_clean_report_ignores_previous_workbook(self):
        from paper_summary import load_inputs
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            from current_release import require_current
            require_current(root,create=True)
            old=root/"le8.paper/dx_comparison.xlsx"
            old.parent.mkdir()
            old.write_bytes(b"old result outside the current analysis")
            data=load_inputs(root,["mi"],["met"],["icd","any1"])
            self.assertTrue(data["availability"].source.eq("no eligible scope").all())
            for key in ("performance","glmnet_neighbors","model_contrasts"):
                self.assertTrue(data[key].empty)

    def test_second_endpoint_runs_on_ordinary_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/"Any1"
            def cmd(first,second):
                return ["endpoint-fits",sys.executable,"-c",f"raise SystemExit({first})","--then",
                        sys.executable,"-c",f"from pathlib import Path;Path({str(out)!r}).touch();raise SystemExit({second})"]
            self.assertEqual(scheduler.main(cmd(9,0)),9)
            self.assertTrue(out.exists())
            out.unlink()
            self.assertEqual(scheduler.main(cmd(0,7)),7)
            self.assertTrue(out.exists())
            out.unlink()
            self.assertEqual(scheduler.main(cmd(143,0)),143)
            self.assertFalse(out.exists())


if __name__=="__main__":unittest.main()
