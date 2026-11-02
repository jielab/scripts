"""Regression tests for actual .cmd allocations and comparison concurrency.

Native index validation and sample discovery are replaced by small fixtures;
config validation, memory budgeting, case generation, .cmd rendering and frozen
plan validation run normally. No scheduler or native classifier is invoked.
"""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpc_cli as cli
import hpc_common as common
import hpc_compare
import hpc_query


class CompareResourcesTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.number = 0
        r1 = self.root / "sample.fastq"
        r1.write_text("@read\nACGT\n+\nIIII\n")
        self.samples = [dict(sample_id="dataset/sample", dataset="dataset", sample="sample",
                             r1=str(r1), r2=None)]
        self.db = self.root / "tiny-index-files"
        self.db.mkdir()
        for name in ("hash.k2d", "taxo.k2d", "opts.k2d", "hash_1.k2d"):
            (self.db / name).write_bytes(b"x" * 64)
        self.taxonomy = self.db / "taxonomy"
        self.taxonomy.mkdir()
        for name in ("nodes.dmp", "names.dmp"):
            (self.taxonomy / name).write_bytes(b"fixture")

    def generate(self, overrides=None, methods=None, repeats=1):
        self.number += 1
        with patch.dict(os.environ, {"MEMORY_DEFAULT_GIB": "160"}):
            cfg = common.load_config()
        cfg.update(work_dir=str(self.root / f"work-{self.number}"), compare_placement="production")
        cfg.update(overrides or {})

        def snapshot(resources):
            engine = resources["classifier_engine"]
            return {"identity": "source" if engine == "kraken" else "converted",
                    "source_database_identity": "source", "db_dir": str(self.db),
                    "files": {name: {} for name in (("hash.k2d", "taxo.k2d", "opts.k2d")
                               if engine == "kraken" else ("hash_1.k2d", "taxo.k2d", "opts.k2d"))},
                    "taxonomy_dir": str(self.taxonomy),
                    "taxonomy_files": {name: {} for name in ("nodes.dmp", "names.dmp")}}

        with patch.object(hpc_query, "discover_samples", return_value=self.samples), \
                patch.object(cli, "engine_snapshot", side_effect=snapshot), \
                contextlib.redirect_stdout(io.StringIO()):
            path = cli.generate_compare(cfg, methods=methods, repeats=repeats, limit=1)
        return cfg, path, json.loads(path.read_text())

    @staticmethod
    def allocation_lines(plan, engine):
        context = plan["context"]
        index = next(i for i, case in enumerate(context["cases"]) if case["classifier_engine"] == engine)
        group = next(name for name, indices in context["groups"].items() if index in indices)
        command = next(c for c in plan["commands"] if c["role"] == "compare_" + group)
        return tuple(line for line in Path(command["path"]).read_text().splitlines()
                     if line.startswith(("#BSUB -q ", "#BSUB -n ", "#BSUB -W ", "#BSUB -R ")))

    def test_independent_memory_and_walltime_reach_submitted_arrays(self):
        _, _, plan = self.generate({"kunpeng_full_memory_gib": 5120, "kunpeng_full_hours": 200},
                                   methods=["kraken", "kunpeng_full"])
        self.assertEqual(len(plan["context"]["groups"]), 2)
        self.assertIn('#BSUB -R "span[hosts=1] rusage[mem=4096GB]"', self.allocation_lines(plan, "kraken"))
        full = self.allocation_lines(plan, "kunpeng_full")
        self.assertIn('#BSUB -R "span[hosts=1] rusage[mem=5120GB]"', full)
        self.assertIn("#BSUB -W 200:00", full)

    def test_independent_queue_is_preserved_and_method_order_does_not_matter(self):
        overrides = {"kunpeng_full_queue": "highmem", "kunpeng_full_memory_gib": 5120}
        _, _, first = self.generate(overrides, methods=["kraken", "kunpeng_full"])
        _, _, reversed_order = self.generate(overrides, methods=["kunpeng_full", "kraken"])
        for engine in ("kraken", "kunpeng_full"):
            self.assertEqual(self.allocation_lines(first, engine), self.allocation_lines(reversed_order, engine))
        self.assertIn("#BSUB -q highmem", self.allocation_lines(first, "kunpeng_full"))
        self.assertIn("#BSUB -q smp", self.allocation_lines(first, "kraken"))

    def test_matching_full_allocations_share_one_group(self):
        _, _, plan = self.generate(methods=["kraken", "kunpeng_full"], repeats=2)
        self.assertEqual(len(plan["context"]["groups"]), 1)
        self.assertEqual(self.allocation_lines(plan, "kraken"), self.allocation_lines(plan, "kunpeng_full"))

    def test_matched_placement_keeps_all_methods_in_one_common_allocation(self):
        _, _, plan = self.generate({"compare_placement": "matched", "kunpeng_full_queue": "highmem",
                                   "kunpeng_full_memory_gib": 5120, "kunpeng_full_hours": 200})
        self.assertEqual(list(plan["context"]["groups"]), ["matched"])
        for engine in ("kunpeng_full", "kunpeng_memory"):
            self.assertEqual(self.allocation_lines(plan, "kraken"), self.allocation_lines(plan, engine))

    def test_split_resource_groups_do_not_multiply_cohort_concurrency(self):
        _, _, plan = self.generate({"kunpeng_full_memory_gib": 5120, "max_array_size": 2}, repeats=3)
        context = plan["context"]
        previous = {}
        all_roles = []
        for group in context["groups"]:
            cohort = context["group_resources"][group]["cohort"]
            role = "compare_" + group
            commands = [c for c in plan["commands"] if c["role"] == role]
            self.assertEqual(len(commands), 2)
            self.assertEqual(commands[0].get("depends_on"), previous.get(cohort))
            if cohort in previous:
                self.assertEqual(commands[0]["dependency_condition"], "ended")
            self.assertEqual(commands[1]["depends_on"], role)
            self.assertEqual(commands[1]["dependency_condition"], "ended")
            previous[cohort] = role
            all_roles.append(role)
        self.assertEqual(set(previous), {"full", "memory"})
        self.assertEqual(plan["commands"][-1]["depends_on"], all_roles)
        self.assertEqual(plan["commands"][-1]["dependency_condition"], "ended")

    def test_worker_rejects_budget_drift_before_starting_native_work(self):
        cfg, path, plan = self.generate(methods=["kraken"])
        group = next(iter(plan["context"]["groups"]))
        resources = common.engine_resources(cfg, "kraken", benchmark=True)
        resources["classify_memory_gib"] += 1
        args = argparse.Namespace(command="compare-worker", allow_local=True, job_plan=str(path),
                                  case_id=1, offset=0, group=group)
        with patch.object(cli, "engine_resources", return_value=resources), \
                patch.object(cli, "engine_snapshot") as snapshot, \
                patch.object(hpc_compare, "run_case") as run_case:
            with self.assertRaisesRegex(RuntimeError, "allocation differs"):
                cli.run_worker(args, cfg)
        snapshot.assert_not_called()
        run_case.assert_not_called()


if __name__ == "__main__":
    unittest.main()
