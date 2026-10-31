import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from abm_core import pipeline
from abm_core.cold_start import preparation_command


class EntryTests(unittest.TestCase):
    def test_abm_stage_runner_preserves_trait_order_and_defaults(self):
        calls = []
        def trait(args, config, name, device):
            calls.append(name)
            self.assertEqual(args.predictor, "blocknet")
            self.assertEqual(args.representation, "csx_ldblocks_fixed")
            self.assertEqual(args.device, "auto")
            self.assertEqual(Path(args.config), pipeline.ROOT / "config/abm.json")
            self.assertEqual(args.through, "report")
        with tempfile.TemporaryDirectory() as analysis_root, \
             patch("abm_core.cold_start.ensure_inputs", side_effect=lambda *a: calls.append("prepare_inputs")), \
             patch("dosage_storage.release_dosages") as release, \
             patch("abm_core.storage.clean_completed_support") as cleanup, \
             patch.object(pipeline, "_backend", return_value=("cpu", None)), \
             patch.object(pipeline, "run_trait", side_effect=trait), \
             patch("final_dispatch.main", side_effect=lambda argv: calls.append("refresh_results")):
            pipeline.main(["--traits", "height,ldl,t2dm", "--analysis-root", analysis_root])
        self.assertEqual(calls, ["prepare_inputs", "height", "ldl", "t2dm"])
        self.assertEqual(cleanup.call_count, 4)
        self.assertEqual(release.call_count, 3)

    def test_dry_run_reveals_defaults_without_touching_data(self):
        text = io.StringIO()
        with patch("abm_core.cold_start.ensure_inputs") as prepare, \
             patch("dosage_storage.release_dosages") as release, \
             patch("abm_core.storage.clean_completed_support") as cleanup, \
             patch.object(pipeline, "run_trait") as trait, \
             patch("final_dispatch.main") as refresh, contextlib.redirect_stdout(text):
            pipeline.main(["--traits", "height,ldl,t2dm", "--dry-run"])
        result = json.loads(text.getvalue())
        self.assertEqual(result["phases"], list(pipeline.PHASES))
        self.assertEqual(result["after_traits"], [])
        self.assertEqual(Path(result["config_path"]).name, "abm.json")
        prepare.assert_not_called(); trait.assert_not_called(); refresh.assert_not_called()
        cleanup.assert_not_called()
        release.assert_not_called()

    def test_partial_stage_does_not_publish_final(self):
        with tempfile.TemporaryDirectory() as analysis_root, \
             patch("abm_core.cold_start.ensure_inputs"), \
             patch("dosage_storage.release_dosages"), \
             patch("abm_core.storage.clean_completed_support"), \
             patch.object(pipeline, "_backend", return_value=("cpu", None)), \
             patch.object(pipeline, "run_trait") as trait, patch("final_dispatch.main") as refresh:
            pipeline.main(["--traits", "height", "--through", "prepare", "--analysis-root", analysis_root])
        self.assertEqual(trait.call_count, 1)
        refresh.assert_not_called()

    def test_missing_inputs_prepare_without_a_source_directory_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = root / "blocks.json"
            config_path.write_text(json.dumps({"prsformer_cache": str(root / "analysis/prsformer/cache")}))
            args = SimpleNamespace(analysis_root=str(root / "analysis"), traits="height", threads=4)
            config = dict(block_config=str(config_path))
            cmd = preparation_command(args, config)
            self.assertEqual(cmd[cmd.index("--prsformer-mode") + 1], "prepare")
            self.assertEqual(cmd[cmd.index("--traits") + 1], "height,ldl,t2dm")
            for flag, relative in (("--cache-dir", "cache/abm"), ("--out-root", "abm"),
                                   ("--prsformer-cache-dir", "prsformer/cache")):
                self.assertEqual(cmd[cmd.index(flag) + 1], str(root / "analysis" / relative))
            self.assertFalse({"--split-file", "--split-group-file", "--keep"} & set(cmd))
            self.assertFalse((root / "analysis").exists())

    def test_prepare_helper_initializes_an_empty_analysis_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            analysis = root / "analysis"
            config_path = root / "blocks.json"
            config_path.write_text(json.dumps({"prsformer_cache": str(analysis / "prsformer/cache")}))
            args = SimpleNamespace(analysis_root=str(analysis), traits="height", threads=4)
            cmd = preparation_command(args, dict(block_config=str(config_path))) + ["--dry-run"]
            env = dict(os.environ, SCRIPT_LOG_DIR=str(root / "logs"), SCRIPT_VERBOSE="1")
            result = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True)
            self.assertIn("main.py --stage prepare", result.stdout)
            self.assertIn(str(analysis / "prsformer/cache/roster"), result.stdout)
            self.assertFalse(analysis.exists())
            split = analysis / "cache/abm/split.rds"
            split.parent.mkdir(parents=True)
            split.write_text("Existing analysis roster; dry-run must not read its contents.")
            result = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True)
            self.assertNotIn("main.py --stage prepare", result.stdout)
            self.assertIn("shared-roster", result.stdout)


if __name__ == "__main__":
    unittest.main()
