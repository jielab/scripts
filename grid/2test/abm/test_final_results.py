import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from abm_core.model import file_sha256, write_json
from abm_core.benchmark_results import performance
from abm_core.prsformer_adapter import publish, verify_complete
import final_compare


TRAITS = ["height", "ldl", "t2dm"]


def fixture(root):
    cache = root / "prsformer/cache"
    cache.mkdir(parents=True)
    data = pd.DataFrame(dict(eid=["001", "002", "003", "004"], split="test",
                             height=[171., 172., np.nan, 174.], ldl=[2., 3., 4., 5.],
                             t2dm=[0., 1., 0., 1.]))
    data.to_csv(cache / "data.tsv.gz", sep="\t", index=False)
    pd.DataFrame(dict(SNP=["rs1"])).to_csv(cache / "variants.tsv.gz", sep="\t", index=False)
    np.save(cache / "genotypes.npy", np.zeros((4, 1), dtype=np.float16))
    meta = dict(signature="fixture", shape=[4, 1],
                sha256={name: file_sha256(cache / name) for name in ("data.tsv.gz", "variants.tsv.gz")})
    write_json(cache / "prepare.json", meta)
    runs = [(trait, run) for trait in TRAITS for run in
            (root / "prsformer" / trait / "run", root / "prsformer2" / trait / "test_run/_LOCAL_ONLY")]
    for trait, run in runs:
        run.mkdir(parents=True)
        (run / "model.pt").write_bytes(b"synthetic model fixture")
        pd.DataFrame(dict(trait=[trait], value=0)).to_csv(run / "metrics.tsv", sep="\t", index=False)
        pd.DataFrame(dict(epoch=[1], loss=[1.])).to_csv(run / "history.tsv", sep="\t", index=False)
        d = data[["eid", "split"]].copy()
        d["trait"], d["y"] = trait, data[trait]
        d["prediction"] = [.2, .7, .3, .8] if trait == "t2dm" else data[trait].fillna(173) + .1
        d.to_csv(run / "test_predictions.tsv.gz", sep="\t", index=False)
        model, genotype = (run / "model.pt").stat(), (cache / "genotypes.npy").stat()
        write_json(run / "evaluation.json", dict(
            format_version=1, preparation_signature="fixture",
            inputs=dict(data_sha256=meta["sha256"]["data.tsv.gz"],
                        variants_sha256=meta["sha256"]["variants.tsv.gz"],
                        genotypes=dict(size=genotype.st_size, mtime_ns=genotype.st_mtime_ns, shape=[4, 1])),
            outputs_sha256={name: file_sha256(run / name) for name in
                            ("metrics.tsv", "history.tsv", "test_predictions.tsv.gz")},
            model=dict(size=model.st_size, mtime_ns=model.st_mtime_ns)))
    for trait in TRAITS:
        run = root / "prsformer2" / trait / "test_run"
        write_json(run / "RUN_MANIFEST.json", {"method": "PRSformer2"})
        (run / "_LOCAL_ONLY/roles.tsv.gz").write_bytes(b"fixture roles")
    return cache


def abm_fixture(root):
    for trait in TRAITS:
        run = root / "abm/FINAL" / trait / "blocknet_fixture"
        run.mkdir(parents=True)
        pd.DataFrame([dict(trait=trait, method="csx-blocknet", n=4)]).to_csv(
            run / "main_performance.tsv", sep="\t", index=False)
        pd.DataFrame([dict(method="csx-blocknet", status="completed")]).to_csv(
            run / "comparator_status.tsv", sep="\t", index=False)
        write_json(run / "RUN_MANIFEST.json", {"contract": {"representation": "csx_ldblocks_fixed"}})
        write_json(run / "COMPLETE.json", dict(status="completed", smoke=False,
                   manifest_sha256=file_sha256(run / "RUN_MANIFEST.json")))
        write_json(run / "PHASES.json", {"phases": {
            "evaluate": {"main_performance.tsv": file_sha256(run / "main_performance.tsv")},
            "report": {"COMPLETE.json": file_sha256(run / "COMPLETE.json")}}})
        write_json(run.parent / "blocknet.CURRENT.json", {"run": run.name})


class FinalTests(unittest.TestCase):
    def test_standalone_final_with_only_abm_and_no_prsformer_inputs(self):
        from final_dispatch import main
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            abm_fixture(root)
            main(["--analysis-root", str(root), "--traits", ",".join(TRAITS)])
            result = pd.read_csv(root / "final/performance.tsv", sep="\t")
            self.assertEqual(set(result.trait), set(TRAITS))
            self.assertEqual(set(result.source), {"ABM"})
            status = pd.read_csv(root / "final/performance_status.tsv", sep="\t")
            missing = status[status.method.isin(["PRSformer", "PRSformer2 exploratory"])]
            self.assertEqual(len(missing), 6)
            self.assertTrue(missing.status.eq("not_completed").all())
            self.assertTrue((root / "final/evidence.sqlite").is_file())
            self.assertFalse((root / "prsformer").exists())
            self.assertFalse((root / "prsformer2").exists())

    def test_final_contains_both_prsformers_without_other_abm_dependencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = fixture(root)
            for trait in TRAITS:
                publish(root / "prsformer2" / trait / "test_run", cache, [trait])
            tables, status = final_compare.current_tables(root, TRAITS)
            result = pd.concat(tables)
            self.assertEqual(set(result.method), {"PRSformer", "PRSformer2 exploratory"})
            self.assertEqual(len(result), 6)
            self.assertEqual(result.loc[result.trait.eq("height"), "n"].tolist(), [3, 3])
            abm_fixture(root)
            final_compare.run(root, TRAITS)
            final = pd.read_csv(root / "final/performance.tsv", sep="\t")
            self.assertEqual(len(final), 9)
            self.assertEqual(set(final.analysis_version), {"FINAL"})
            self.assertEqual(set(final.method), {"PRSformer", "PRSformer2 exploratory", "csx-blocknet"})
            self.assertTrue((root / "final/performance.xlsx").is_file())
            from browser_index import Index
            result = Index(root, root / "final").run(TRAITS)
            self.assertGreater(result["tables"], 0)
            catalogue = pd.read_csv(root / "final/catalogue.tsv", sep="\t")
            self.assertIn("prsformer2", set(catalogue.module))
            self.assertTrue(catalogue.analysis_version.eq("FINAL").all())
            self.assertFalse(catalogue.source.str.contains("_LOCAL_ONLY").any())

    def test_final_rejects_changed_prsformer_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = fixture(root)
            for trait in TRAITS:
                publish(root / "prsformer2" / trait / "test_run", cache, [trait])
            (root / "prsformer/height/run/test_predictions.tsv.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "changed after completion"):
                final_compare.current_tables(root, TRAITS)

    def test_prsformer2_receipt_and_abm_reports_reject_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = fixture(root)
            run = root / "prsformer2/height/test_run"
            publish(run, cache, ["height"])
            verify_complete(run, cache, ["height"])
            (run / "main_performance.tsv").write_text("changed")
            with self.assertRaisesRegex(ValueError, "changed"):
                verify_complete(run, cache, ["height"])
            abm_fixture(root)
            (root / "abm/FINAL/height/blocknet_fixture/main_performance.tsv").write_text("changed")
            with self.assertRaisesRegex(ValueError, "report changed"):
                final_compare.current_tables(root, TRAITS)

    def test_test_roster_and_outcomes_are_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = fixture(root)
            # Accept the current prepared cohort without a hard-coded historical hash.
            performance(root / "prsformer/height/run", cache, ["height"], "PRSformer")
            path = root / "prsformer/height/run/test_predictions.tsv.gz"
            d = pd.read_csv(path, sep="\t", dtype={"eid": str})
            receipt_path = path.parent / "evaluation.json"
            receipt = json.loads(receipt_path.read_text())
            for column, value, message in (("eid", "not-in-current-cohort", "incomplete test roster"),
                                            ("y", d.loc[0, "y"] + 1, "outcomes differ")):
                with self.subTest(column=column):
                    changed = d.copy()
                    changed.loc[0, column] = value
                    changed.to_csv(path, sep="\t", index=False)
                    receipt["outputs_sha256"][path.name] = file_sha256(path)
                    write_json(receipt_path, receipt)
                    with self.assertRaisesRegex(ValueError, message):
                        performance(path.parent, cache, ["height"], "PRSformer")
