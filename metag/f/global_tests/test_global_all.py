"""Default batch orchestration and failure boundaries, without production data."""
from contextlib import ExitStack, redirect_stdout, redirect_stderr
import io
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

import global_pipeline as cli
import global_build
import global_query


class AllStagesTests(unittest.TestCase):
    def test_default_explicit_all_and_help(self):
        self.assertEqual(cli.parse_args([]).command, "all")
        for prefix in ([], ["all"]):
            args = cli.parse_args(prefix + ["--dat-list", "changj,Ndesert,Wxiamen",
                                           "--threads", "3", "--benchmark-fragments", "17"])
            self.assertEqual(args.command, "all")
            self.assertEqual(args.dat_list, "changj,Ndesert,Wxiamen")
            self.assertEqual(args.benchmark_fragments, 17)
            self.assertEqual(args.db_gib, 24)
            self.assertEqual(args.memory_limit_gib, 40)
            self.assertTrue(args.quick)
            self.assertFalse(args.read_reports)
            self.assertFalse(args.no_masking)
        with redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as result:
            cli.main(["--help"])
        self.assertEqual(result.exception.code, 0)
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            cli.parse_args(["splti"])

    def exercise(self, fail=None, custom=False):
        events = []
        samples = [{"name": "changj__first"}, {"name": "Ndesert__second"}]
        def record(name, value=None):
            def call(*args, **kwargs):
                events.append(name)
                if fail == name:
                    raise RuntimeError(name)
                if name == "build":
                    self.assertEqual(args[0].chunk_bases, 3 << 20)
                    self.assertEqual(args[0].threads, 7)
                if name == "benchmark":
                    self.assertEqual(args[0], samples[0])
                    self.assertEqual(args[3].benchmark_fragments, 17)
                    self.assertEqual(args[3].confidence, 0.2)
                return value
            return call
        with ExitStack() as stack:
            stack.enter_context(redirect_stdout(io.StringIO()))
            stack.enter_context(patch.dict(cli.os.environ, {"METAG_KRAKEN2_DIR": ""}))
            for module, name, action in (
                (cli, "assert_reference_policy", lambda: None),
                (cli, "ensure_memory_scope", lambda *args: None),
                (cli.subprocess, "run", record("install")),
                (cli, "configure_tools", record("tools", Path("/fixture/kraken2"))),
                (cli, "plan", record("plan", 2 if fail == "plan_status" else 0)),
                (cli, "discover_samples", lambda args: samples),
                (global_build, "build_database", record("build")),
                (global_build, "database_snapshot", record("snapshot", {"path": "fixture"})),
                (global_query, "benchmark", record("benchmark")),
                (global_query, "run_sample", lambda sample, *args: events.append(sample["name"])),
            ):
                stack.enter_context(patch.object(module, name, side_effect=action))
            argv = ["--dat-list", "changj,Ndesert", "--chunk-mib", "3",
                    "--threads", "7", "--benchmark-fragments", "17", "--confidence", "0.2"]
            if custom:
                argv += ["--kraken-dir", "/fixture"]
            if fail and fail != "plan_status":
                with self.assertRaisesRegex(RuntimeError, fail):
                    cli.main(argv)
            else:
                self.assertEqual(cli.main(argv), 2 if fail else 0)
        return events

    def test_default_runs_all_stages_and_all_selected_samples_in_order(self):
        self.assertEqual(self.exercise(), ["install", "plan", "tools", "build", "snapshot",
                                          "benchmark", "changj__first", "Ndesert__second"])

    def test_every_failure_stops_downstream_work(self):
        expected = ["install", "plan", "tools", "build", "snapshot", "benchmark"]
        for index, stage in enumerate(expected):
            with self.subTest(stage=stage):
                self.assertEqual(self.exercise(fail=stage), expected[:index + 1])
        self.assertEqual(self.exercise(fail="plan_status"), ["install", "plan"])

    def test_explicit_installation_is_verified_without_default_install(self):
        events = self.exercise(custom=True)
        self.assertNotIn("install", events)
        self.assertEqual(events[0], "tools")


if __name__ == "__main__":
    unittest.main()
