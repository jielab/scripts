"""Meaningful input-integrity checks; native CIGMA estimation is NOT exercised."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"f"))
from c3_cigma import read_inputs, bh, matrix
from c3_cell_annotation import annotate


class InputIntegrity(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.rng = np.random.default_rng(2026)
        self.ids = [f"{i:03d}" for i in range(24)]
        self.cells = ["monocyte","hepatocyte"]
        x = self.rng.normal(size=(24,6))
        self.k = x@x.T/6 + np.eye(24)*.2
        self.y = self.rng.normal(size=(24,2))
        for name, values, cols in (("ctp",self.y,self.cells),("ctnu",np.ones((24,2))*.01,self.cells),
                                  ("P",np.ones((24,2))*.5,self.cells),("K",self.k,self.ids)):
            pd.DataFrame(values,index=self.ids,columns=cols).rename_axis("eid").to_csv(self.root/f"{name}.csv")
        self.row = dict(ctp="ctp.csv",ctnu="ctnu.csv",P="P.csv",K="K.csv",ctnu_definition="variance_of_pseudobulk_mean")

    def tearDown(self): self.tmp.cleanup()

    def test_alignment_uses_labels_and_preserves_leading_zero_ids(self):
        for name in ("ctnu","P","K"):
            d = matrix(self.root/f"{name}.csv").iloc[::-1,::-1]
            d.to_csv(self.root/f"{name}.csv")
        args,cells,_ = read_inputs(self.row,self.root)
        np.testing.assert_allclose(args["K"],self.k)
        np.testing.assert_allclose(args["Y"],self.y)
        self.assertEqual(cells,self.cells)

    def test_wrong_noise_units_are_rejected(self):
        for unit in ("SD","cell_variance",""):
            with self.assertRaises(ValueError):read_inputs(dict(self.row,ctnu_definition=unit),self.root)

    def test_unidentifiable_identity_kinship_rejected(self):
        pd.DataFrame(np.eye(24),index=self.ids,columns=self.ids).to_csv(self.root/"K.csv")
        with self.assertRaises(ValueError):read_inputs(self.row,self.root)

    def test_missing_donor_rejected(self):
        d = matrix(self.root/"P.csv");d.index = ["different"]+self.ids[1:];d.to_csv(self.root/"P.csv")
        with self.assertRaises(ValueError):read_inputs(self.row,self.root)

    def test_missing_hypotheses_retained_in_FDR_family(self):
        np.testing.assert_allclose(bh([.01,.04,np.nan])[:2],[.03,.06])
        with self.assertRaises(ValueError):bh([1.2])

    def test_two_assays_same_gene_not_two_enrichment_hits(self):
        pd.DataFrame(dict(assay=["NPPB","NTPROBNP","PCSK9","IL6"],gene=["NPPB","NPPB","PCSK9","IL6"])).to_csv(self.root/"u.csv",index=False)
        pd.DataFrame(dict(gene=["NPPB","PCSK9"],cell_type=["cardiomyocyte","hepatocyte"],source=["synthetic fixture"]*2,source_version=["test"]*2)).to_csv(self.root/"a.csv",index=False)
        pd.DataFrame(dict(model=["panel"]*2,feature=["NPPB","NTPROBNP"])).to_csv(self.root/"p.csv",index=False)
        annotate(self.root/"u.csv",self.root/"a.csv",self.root/"p.csv",self.root/"out")
        z = pd.read_csv(self.root/"out/cell.enrichment.csv")
        self.assertTrue((z.background_genes==3).all())
        self.assertTrue((z.panel_genes==1).all())
        self.assertEqual(int(z.loc[z.cell_type=="cardiomyocyte","hits"].iloc[0]),1)


if __name__ == "__main__": unittest.main()
