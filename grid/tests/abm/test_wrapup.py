import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd

from abm_core import cold_start, prsformer_adapter
import final_dispatch
from final_outputs import PAPER_FILES, aggregate
import paper_summary
import prsformer_prepare_shared
from prepare_kinship import components

ROOT = Path(__file__).resolve().parents[2]


class WrapupTests(unittest.TestCase):
    def shell(self, script, *args):
        with tempfile.TemporaryDirectory() as tmp:
            return subprocess.run([str(ROOT / script), *args], text=True, capture_output=True,
                                  env=dict(os.environ, GRID_PYTHON=sys.executable,
                                           SCRIPT_LOG_DIR=tmp, SCRIPT_VERBOSE="1"), check=True).stdout

    def test_public_module_lists_and_default_benchmarks(self):
        text = self.shell("grid.sh", "--help")
        modules = text.split("Examples:")[0].lower()
        for module in ("blocknet", "twins", "assess", "final"):
            self.assertIn(module, modules)
        for module in ("prsformer", "prsformer2", "workflow"):
            self.assertNotIn(module, modules)
        text = self.shell("4.prsformer.sh", "--help")
        for module in ("prepare", "prsformer", "prsformer2", "all"):
            self.assertIn(module, text)
        text = self.shell("4.prsformer.sh", "--traits", "height,ldl,t2dm", "--dry-run")
        first, second, third = (text.index(s) for s in
                                ("[4.prsformer.sh prepare]", "PRSformer plan:", "PRSformer2 plan:"))
        self.assertLess(first, second)
        self.assertLess(second, third)

    def test_named_phases_are_existing_pipeline_views(self):
        for module, phase in (("blocknet", "fit_predictors"), ("twins", "build_support"),
                              ("assess", "report"), ("abm", "report")):
            with self.subTest(module=module):
                text = self.shell("grid.sh", module, "--traits", "height", "--dry-run")
                plan, _ = json.JSONDecoder().raw_decode(text[text.index('{'):])
                self.assertEqual(plan["phases"][-1], phase)
        text = self.shell("grid.sh", "--traits", "height", "--dry-run")
        self.assertIn("paper_summary.py", text)
        self.assertNotIn("PRSformer plan:", text)

    def test_prepare_and_final_dry_runs_have_no_side_effects(self):
        with tempfile.TemporaryDirectory() as tmp, patch("subprocess.run") as run, \
             patch("abm_prepare.prepare_abm") as prepare, contextlib.redirect_stdout(io.StringIO()):
            root = str(Path(tmp) / "absent")
            prsformer_prepare_shared.main(["--analysis-root", root, "--dry-run"])
            final_dispatch.main(["--analysis-root", root, "--dry-run"])
            run.assert_not_called()
            prepare.assert_not_called()
            self.assertFalse(Path(root).exists())

    def test_missing_inputs_never_start_preparation(self):
        with tempfile.TemporaryDirectory() as tmp, patch("subprocess.run") as run:
            root = Path(tmp)
            cfg = root / "blocks.json"
            cfg.write_text(json.dumps(dict(prsformer_cache=str(root / "prsformer/cache"))))
            args = SimpleNamespace(analysis_root=tmp, traits="height", prepared=None)
            with self.assertRaisesRegex(FileNotFoundError, "4.prsformer.sh prepare"):
                cold_start.ensure_inputs(args, dict(block_config=str(cfg)))
            run.assert_not_called()

    def test_bound_blocks_are_verified_without_genotypes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = root / "blocks.json"
            config = dict(prsformer_cache=str(root / "cache"), ld_block_dir=tmp, ld_block_files="x")
            cfg.write_text(json.dumps(config))
            marker = root / "BLOCKS_COMPLETE.json"
            marker.touch()
            args = SimpleNamespace(analysis_root=tmp, traits="height", predictor="blocknet")
            array = MagicMock()
            with patch("block_inputs.read_config", return_value=config), \
                 patch("block_inputs.resolve_block_files", return_value={}), \
                 patch("block_inputs.cache_paths", return_value=[marker]), \
                 patch("abm_core.grid_inputs.locate_blocks", return_value=[array]) as verify:
                cold_start.ensure_inputs(args, dict(block_config=str(cfg)))
                self.assertFalse(verify.call_args.kwargs["allow_build"])
                array._mmap.close.assert_called_once()
                verify.side_effect = ValueError("untrusted producer")
                with self.assertRaisesRegex(ValueError, "untrusted producer"):
                    cold_start.ensure_inputs(args, dict(block_config=str(cfg)))

    def test_prsformer2_wrapper_environment_retains_dosage(self):
        with patch.dict(os.environ, GRID_KEEP_SHARED_INPUTS="1"), \
             patch.object(prsformer_adapter, "_run") as run, \
             patch.object(prsformer_adapter, "benchmark_device", return_value="cuda"), \
             patch("dosage_storage.release_dosages") as release:
            prsformer_adapter.main(["--traits", "height"])
            run.assert_called_once()
            release.assert_not_called()

    def test_final_dispatch_order_and_no_cleanup(self):
        with patch("subprocess.run") as run, patch("cache_cleanup.clean_csx") as clean:
            final_dispatch.main([])
        self.assertEqual([Path(call.args[0][1]).name for call in run.call_args_list],
                         ["final_compare.py", "paper_summary.py", "browser_index.py"])
        clean.assert_not_called()

    def test_empty_or_incomplete_benchmarks_replace_stale_paper_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "final"
            out.mkdir()
            (out / "paper_main_performance.tsv").write_text("trait\tmethod\nheight\tstale\n")
            parent = root / "prsformer2/height"
            parent.mkdir(parents=True)
            (parent / "CURRENT.json").write_text('{"run": "unfinished"}')
            (parent / "unfinished").mkdir()
            with contextlib.redirect_stdout(io.StringIO()):
                paper_summary.main(["--analysis-root", tmp, "--traits", "height"])
            for name in PAPER_FILES:
                aggregate(pd.read_csv(out / name, sep="\t"))
            self.assertEqual(len(pd.read_csv(out / "paper_main_performance.tsv", sep="\t")), 0)
            status = pd.read_csv(out / "paper_method_status.tsv", sep="\t")
            self.assertEqual(status.loc[status.component.eq("PRSformer2 exploratory"), "status"].item(), "not_completed")

    def test_identifier_columns_never_enter_paper_outputs(self):
        for name in ("eid", "IID", "participant_id", "donor.id", "group_id", "prediction", "y"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Individual-level"):
                aggregate(pd.DataFrame({name: ["private"]}))
        aggregate(pd.DataFrame(dict(run_id=["run"], prediction_r2=[.2], n_groups=[20])))

    def test_paper_uses_exact_current_and_rejects_tampering(self):
        from test_final_results import abm_fixture
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            abm_fixture(root)
            parent = root / "abm/FINAL/height"
            pointer = parent / "blocknet.CURRENT.json"
            pointer.rename(parent / "blocknet_old.CURRENT.json")
            self.assertIsNone(paper_summary.current_blocknet(root, "height"))
            (parent / "blocknet_old.CURRENT.json").rename(pointer)
            run = paper_summary.current_blocknet(root, "height")
            (run / "main_performance.tsv").write_text("changed")
            with self.assertRaisesRegex(ValueError, "report changed"):
                paper_summary.current_blocknet(root, "height")

    def test_components_keep_outside_cohort_paths(self):
        edges = pd.DataFrame(dict(ID1=["a", "outside"], ID2=["outside", "b"]))
        got = components(edges, ["a", "b", "singleton"])
        self.assertEqual(got.iloc[0], got.iloc[1])
        self.assertNotEqual(got.iloc[0], got.iloc[2])


if __name__ == "__main__":
    unittest.main()
