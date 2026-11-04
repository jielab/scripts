"""HPC control-plane tests. No real LSF jobs or scientific data are submitted."""
import argparse
import contextlib
import csv
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("gu_hpc", Path(__file__).resolve().parents[2] / "f" / "hpc.py")
hpc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hpc)


class HPC(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gu-hpc-test-")
        self.base = Path(self.tmp.name)
        self.root = self.base / "code with spaces" / "gu"
        self.root.mkdir(parents=True)
        (self.root / "gu.sh").write_text("#!/bin/bash\nexit 0\n")
        (self.root / "arg.sh").write_text("#!/bin/bash\nexit 0\n")
        self.saved_root = hpc.ROOT
        hpc.ROOT = self.root
        self.env = self.base / "runtime.env"
        self.env.write_text("export EXAMPLE=1\n")
        self.gen = self.base / "pfile"
        self.gen.mkdir()
        for c in ["1", "2", "X"]:
            for ext, content in [("pgen", "genotype placeholder"), ("pvar", "#CHROM\tPOS\tID\tREF\tALT\n"),
                                 ("psam", "#IID SEX\na 1\nb 2\nc 1\nd 2\ne 1\n")]:
                (self.gen / f"chr{c}.{ext}").write_text(content)
        self.run = self.base / "run with 'quote'"
        # Tests use private scratch; actual CLI refuses ephemeral result roots.
        self.persistent = patch.object(hpc, "persistent_path", side_effect=lambda p: Path(p).resolve())
        self.persistent.start()

    def tearDown(self):
        self.persistent.stop()
        hpc.ROOT = self.saved_root
        self.tmp.cleanup()

    def make_plan(self, method="ibdmix", *more):
        values = ["plan", method, "--run-dir", str(self.run), "--env-file", str(self.env),
                  "--queue", "ser", "--dir-target", str(self.gen / "chr"), "--chr", "1,2,X"]
        if method in {"ibdmix", "phyml"}:
            values += ["--archaic-path", str(self.base / "archaic")]
        if method in {"arg", "arg-trace"}:
            values += ["--dir-gen", str(self.gen.parent)]
        if method == "trace":
            values += ["--dir-arg", str(self.base / "trees")]
        args = hpc.parser().parse_args(values + list(more))
        with contextlib.redirect_stdout(io.StringIO()):
            hpc.plan_run(args)
        return hpc.load_plan(self.run)

    def output(self, task, rows=1, archive=False):
        p = Path(task["output"]) / "ibdmix" / ("chr" + task["chromosomes"][0]) / "final" / "all_archaic_refs.lod4.len50000.segments.tsv.gz"
        p.parent.mkdir(parents=True)
        with gzip.open(p, "wt") as f:
            f.write("ID\tchrom\tstart\tend\tanc\n")
            for n in range(rows):
                f.write(f"s{n}\t{task['chromosomes'][0]}\t100\t50100\tAltai\n")
        if archive:
            bundle = p.parents[1] / "ibdmix.raw.tar.gz"
            with tarfile.open(bundle, "w:gz") as t:
                t.add(p, arcname="final/" + p.name)
        artifacts = hpc.artifact_candidates(task)
        hpc.validate_artifacts(task, artifacts)
        return artifacts

    def seal(self, plan, task, artifacts):
        hpc.atomic_json(hpc.task_receipt(plan, task), {"plan_id": plan["plan_id"], "state": "complete", "artifacts": artifacts})

    def test_chromosome_ranges(self):
        self.assertEqual(hpc.chroms("chr1,2-4,23,chrX,1"), ["1", "2", "3", "4", "X"])
        for invalid in ["0", "2-1", "Y", "1-23", ""]:
            with self.assertRaises(ValueError):
                hpc.chroms(invalid)

    def test_batch_then_sex_filter(self):
        rows = hpc.read_psam(self.gen / "chr1.psam")
        batches = hpc.freeze_batches(rows, 2, True)
        self.assertEqual([b["scored_n"] for b in batches], [1, 1, 1])
        self.assertEqual([b["rows"][0][1] for b in batches], ["a", "c", "e"])
        self.assertEqual([b["first_row"] for b in batches], [1, 3, 5])

    def test_duplicate_iid_rejected(self):
        p = self.base / "duplicate.psam"
        p.write_text("#FID IID SEX\n0 a 1\n1 a 2\n")
        with self.assertRaisesRegex(ValueError, "Duplicate IID"):
            hpc.read_psam(p)

    def test_batch_chromosome_isolation_and_cpu(self):
        plan = self.make_plan("ibdmix", "--psam", str(self.gen / "chr1.psam"), "--batch-size", "2", "--keep-males")
        self.assertEqual(len(plan["tasks"]), 9)
        self.assertEqual(len({t["output"] for t in plan["tasks"]}), 9)
        for task in plan["tasks"]:
            args = task["command"]
            self.assertEqual(args[args.index("--GPU") + 1], "FALSE")
            self.assertEqual(args[args.index("--jobs") + 1], "1")
            self.assertEqual(args[args.index("--foreground") + 1], "TRUE")
            self.assertNotIn("--keep-psam", args)
        self.assertEqual(Path(plan["batches"][0]["keep"]).read_text(), "#FID\tIID\n0\ta\n")

    def test_trace_keeps_joint_chromosomes(self):
        plan = self.make_plan("trace")
        self.assertEqual(len(plan["tasks"]), 1)
        self.assertEqual(plan["tasks"][0]["chromosomes"], ["1", "2", "X"])
        cmd = plan["tasks"][0]["command"]
        self.assertEqual(cmd[cmd.index("--chr") + 1], "1,2,X")

    def test_trace_batch_refused(self):
        with self.assertRaisesRegex(ValueError, "fixed complete cohort"):
            self.make_plan("trace", "--psam", str(self.gen / "chr1.psam"))

    def test_arg_trace_dependency_plan(self):
        plan = self.make_plan("arg-trace", "--memory-gb", "96")
        self.assertEqual(len(plan["tasks"]), 3)
        self.assertEqual(plan["trace_task"]["kind"], "trace")
        for task in plan["tasks"]:
            cmd = task["command"]
            self.assertEqual(cmd[cmd.index("--memory-limit-gb") + 1], "96")
            self.assertIn("--dir-gen4arg", cmd)
            self.assertIn("--dir-arg", cmd)
            self.assertEqual(cmd[cmd.index("--dir-map") + 1], "/work/sph-huangj/data/gen/hm3/37/gmap")

    def test_fresh_reference_defaults_and_cli(self):
        plan = self.make_plan()
        self.assertEqual(plan["archaic_gen"], "Altai,Chagyr,Denisova,Denisova25,Vindija")
        values = ["plan", "trace", "--run-dir", str(self.run), "--env-file", str(self.env),
                  "--queue", "ser", "--dir-target", str(self.gen / "chr")]
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            hpc.parser().parse_args(values + ["--arg-dir", "/data/trees"])

    def test_arg_cohort_uses_frozen_sample_panel(self):
        panel = self.base / "samples.txt"
        panel.write_text("sample\tpop\tsuper_pop\tsex\na\tCEU\tEUR\t1\n")
        plan = self.make_plan("arg-trace", "--sample-panel", str(panel))
        for task in plan["tasks"]:
            cmd = task["command"]
            self.assertEqual(cmd[cmd.index("--ancestry-file") + 1], plan["sample_panel"])
            self.assertEqual(cmd[cmd.index("--sample-panel") + 1], plan["sample_panel"])

    def test_runtime_configuration_is_part_of_plan_identity(self):
        runtime = self.root / "f" / "runtime.env"
        runtime.parent.mkdir()
        runtime.write_text('export GU_DATA_ROOT="/work/sph-huangj/data"\n')
        plan = self.make_plan()
        self.assertIn(str(runtime), plan["source_sha256"])
        runtime.write_text('export GU_DATA_ROOT="/different/data"\n')
        with self.assertRaisesRegex(ValueError, "Code changed"):
            hpc.verify_plan_files(plan)

    def test_deployed_commands_use_hpc_paths(self):
        plan = self.make_plan()
        plan.update(source_root="/work/sph-huangj/scripts.HPC/gu",
                    dir_target="/work/sph-huangj/data/gen/ukb/37/imp/chr", target="ukb",
                    archaic_path="/work/sph-huangj/data/gen/archaic/37/vcf",
                    environment="/work/sph-huangj/data/analysis/gu/test/inputs/runtime.env")
        with patch.object(hpc, "ROOT", Path(plan["source_root"])):
            command = hpc.common_gu(plan, "/work/sph-huangj/data/analysis/gu/test/tasks/batch00001/chr1")
            plan["tasks"][0]["command"] = command
            hpc.render_jobs(plan)
        self.assertEqual(command[1], "/work/sph-huangj/scripts.HPC/gu/gu.sh")
        rendered = "\n".join(p.read_text() for p in (self.run / "jobs").glob("*.cmd"))
        self.assertIn("/work/sph-huangj/scripts.HPC/gu/f/hpc.py", rendered)
        for old_mount in ("/mnt/f", "/mnt/d", "/workspace/scratch"):
            self.assertNotIn(old_mount, rendered + "\n" + " ".join(command))

    def test_plan_tamper_rejected(self):
        plan = self.make_plan()
        plan["threads"] = 999
        hpc.atomic_json(self.run / "plan.json", plan)
        with self.assertRaisesRegex(ValueError, "edited"):
            hpc.load_plan(self.run)

    def test_input_and_code_drift_rejected(self):
        plan = self.make_plan()
        hpc.verify_plan_files(plan)
        (self.gen / "chr1.psam").write_text("#IID SEX\nz 1\n")
        with self.assertRaisesRegex(ValueError, "Input changed"):
            hpc.verify_plan_files(plan)

    def test_shell_quoting_and_login_guard(self):
        self.make_plan()
        for p in (self.run / "jobs").glob("*.cmd"):
            result = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        env = dict(os.environ)
        env.pop("LSB_JOBID", None)
        result = subprocess.run(["bash", str(self.run / "jobs" / "worker.cmd")], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("login node", result.stderr)

    def test_scheduler_slots_guard(self):
        plan = self.make_plan()
        with patch.dict(os.environ, {"LSB_JOBID": "123", "LSB_DJOB_NUMPROC": "1"}):
            with self.assertRaisesRegex(ValueError, "slots"):
                hpc.require_scheduler(plan, "worker")
            hpc.require_scheduler(plan, "prepare")

    def test_short_queue_minimum_for_every_stage(self):
        plan = self.make_plan()
        plan['queue'] = 'short'
        for stage in ('prepare', 'worker', 'merge', 'trace'):
            cmd = hpc.bsub_args(plan, stage, 'gu.test', None)
            self.assertEqual(cmd[cmd.index('-n') + 1], '40')
        plan['slots'] = 65
        cmd = hpc.bsub_args(plan, 'worker', 'gu.test', None)
        self.assertEqual(cmd[cmd.index('-n') + 1], '65')

    def test_qiming_budget_counts_the_supervisor_slot(self):
        with self.assertRaisesRegex(ValueError, "121 cores"):
            self.make_plan("ibdmix", "--queue", "1t81c", "--threads", "120", "--memory-gb", "768")
        plan = self.make_plan("ibdmix", "--queue", "1t81c", "--threads", "119", "--memory-gb", "768")
        self.assertEqual(plan["slots"], 120)

    def test_flat_manifest_includes_helpers_excludes_other_projects(self):
        (self.root / 'f').mkdir()
        helper = self.root / 'f/results.py'
        helper.write_text('original')
        (self.root / 'ld').mkdir()
        (self.root / 'ld/ld.sh').write_text('unrelated')
        before = hpc.code_manifest()
        self.assertIn(str(helper), before)
        self.assertNotIn(str(self.root / 'ld/ld.sh'), before)
        helper.write_text('changed')
        self.assertNotEqual(before[str(helper)], hpc.code_manifest()[str(helper)])

    def test_memory_units_and_array_name(self):
        plan = self.make_plan()
        cmd = hpc.bsub_args(plan, "worker", "gu.unique.a", "ended(123)", 80, 10)
        self.assertEqual(cmd[cmd.index("-M") + 1], "32GB")
        self.assertIn("rusage[mem=32768MB/host]", cmd[cmd.index("-R") + 1])
        self.assertEqual(cmd[cmd.index("-J") + 1], "gu.unique.a[1-80]%10")
        self.assertEqual(cmd[cmd.index("-w") + 1], "ended(123)")
        plan["rusage_memory_mode"] = "per-slot"
        cmd = hpc.bsub_args(plan, "worker", "x", None)
        self.assertIn("rusage[mem=6554MB/task]", cmd[cmd.index("-R") + 1])
        plan["rusage_memory_mode"] = "per-job"
        cmd = hpc.bsub_args(plan, "worker", "x", None)
        self.assertIn("rusage[mem=32768MB/job]", cmd[cmd.index("-R") + 1])

    def test_gpu_backend_and_resource_request_are_paired(self):
        plan = self.make_plan("ibdmix", "--GPU", "TRUE", "--queue", "4a100-40")
        self.assertTrue(plan["gpu"])
        for task in plan["tasks"]:
            command = task["command"]
            self.assertEqual(command[command.index("--GPU") + 1], "TRUE")
        worker = hpc.bsub_args(plan, "worker", "gpu", "done(100)")
        self.assertEqual(worker[worker.index("-q") + 1], "4a100-40")
        self.assertIn("num=1/host:mode=exclusive_process:j_exclusive=yes", worker)
        for stage in ("prepare", "merge"):
            command = hpc.bsub_args(plan, stage, "cpu", None)
            self.assertNotIn("-gpu", command)
            self.assertEqual(command[command.index("-q") + 1], plan["cpu_queue"])
        self.assertEqual(hpc.runtime_environment(plan)["IBDMIX_BACKEND"], "cuda")

    def test_phyml_timeout_policy_is_frozen_in_plan(self):
        loci = self.base / "leads.jma.cojo"
        loci.write_text("Chr SNP bp\n1 rs1 1\n")
        plan = self.make_plan("phyml", "--loci", str(loci), "--phyml-timeout-seconds", "86400", "--phyml-retry-timeouts")
        env = hpc.runtime_environment(plan)
        self.assertEqual(env["PHYML_TREE_TIMEOUT"], "86400")
        self.assertEqual(env["PHYML_RETRY_TIMEOUTS"], "TRUE")

    def test_saved_plan_freezes_super_roots_and_shared_helpers(self):
        shared = self.root.parent / "0f"
        shared.mkdir()
        os_config = shared / "0os.sh"
        os_config.write_text("# original path policy\n")
        overrides = {"DIR_DAT":"/work/custom small", "DIR_BIG":"/work/custom-big", "DIR_TMP":"/scratch/custom",
                     "DIR_ANNOT":"/work/custom-annotation", "GU_SOFT":"/work/custom-software"}
        with patch.dict(os.environ, overrides):
            plan = self.make_plan()
        self.assertTrue(overrides.items() <= plan["resolved_environment"].items())
        self.assertIn(str(os_config), plan["source_sha256"])
        text = (self.run / "jobs" / "worker.cmd").read_text()
        self.assertIn("export DIR_DAT='/work/custom small'", text)
        self.assertLess(text.index("export DIR_TMP="), text.index("source "))
        os_config.write_text("# changed path policy\n")
        with self.assertRaisesRegex(ValueError, "Code changed"):
            hpc.verify_plan_files(plan)

    def test_archive_and_materialized_not_double_counted(self):
        plan = self.make_plan()
        task = plan["tasks"][0]
        artifacts = self.output(task, rows=2, archive=True)
        self.assertEqual(len(artifacts), 1)
        self.assertEqual(artifacts[0]["rows"], 2)
        self.assertIn("member", artifacts[0])

    def test_empty_success_distinct_from_missing(self):
        plan = self.make_plan()
        for task in plan["tasks"]:
            self.seal(plan, task, self.output(task, rows=0))
        hpc.merge(plan)
        rec = hpc.read_json(self.run / "summary" / "COMPLETE.json")
        self.assertEqual(rec["rows"], 0)
        with gzip.open(self.run / "summary" / "segments.tsv.gz", "rt") as f:
            self.assertEqual(len(f.readlines()), 1)

    def test_partial_results_block_merge(self):
        plan = self.make_plan()
        self.seal(plan, plan["tasks"][0], self.output(plan["tasks"][0]))
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            hpc.merge(plan)
        self.assertFalse((self.run / "summary" / "COMPLETE.json").exists())

    def test_changed_results_block_merge(self):
        plan = self.make_plan()
        for task in plan["tasks"]:
            self.seal(plan, task, self.output(task))
        file = next(Path(plan["tasks"][0]["output"]).rglob("*.segments.tsv.gz"))
        file.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "outputs_changed"):
            hpc.merge(plan)

    def test_aggregate_retains_batch_provenance(self):
        plan = self.make_plan()
        for i, task in enumerate(plan["tasks"]):
            self.seal(plan, task, self.output(task, rows=i))
        hpc.merge(plan)
        with gzip.open(self.run / "summary" / "segments.tsv.gz", "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        self.assertEqual(len(rows), 3)
        self.assertIn("hpc_batch", rows[0])
        self.assertNotIn("ancestry_percent", rows[0])

    def test_mock_lsf_submission_and_duplicate_guard(self):
        plan = self.make_plan("ibdmix", "--array-size", "2", "--max-running", "3")
        calls = []
        def mock_run(cmd, **kwargs):
            calls.append((cmd, kwargs))
            return subprocess.CompletedProcess(cmd, 0, f"Job <{700 + len(calls)}> is submitted to queue <ser>.\n", "")
        with patch.object(hpc.shutil, "which", return_value="/mock/bsub"), patch.object(hpc.subprocess, "run", side_effect=mock_run):
            with contextlib.redirect_stdout(io.StringIO()):
                hpc.submit(argparse.Namespace(run_dir=str(self.run)))
            self.assertEqual(len(calls), 4)
            self.assertEqual(calls[1][0][calls[1][0].index("-w") + 1], "ended(701)")
            self.assertEqual(calls[-1][0][calls[-1][0].index("-w") + 1], "ended(702) && ended(703)")
            self.assertIn("[1-2]%2", calls[1][0][calls[1][0].index("-J") + 1])
            with self.assertRaisesRegex(ValueError, "already submitted"):
                hpc.submit(argparse.Namespace(run_dir=str(self.run)))

    def test_active_retry_blocked(self):
        events = [{"state": "submitted", "stage": "worker", "job_id": "123", "selection": [1]}]
        with patch.object(hpc.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "RUN\n", "")):
            with self.assertRaisesRegex(ValueError, "still active"):
                hpc.assert_inactive(events)

    def test_batch_id_is_exported_to_scientific_worker(self):
        plan = self.make_plan("ibdmix", "--psam", str(self.gen / "chr1.psam"), "--batch-size", "2")
        task = next(t for t in plan["tasks"] if t["batch"] == 2)
        hpc.atomic_json(self.run / "state" / "prepare.json", {"plan_id": plan["plan_id"], "state": "complete"})
        environments = []
        def scientific_stub(command, **kwargs):
            environments.append(kwargs["env"])
            self.output(task, rows=0)
            return subprocess.CompletedProcess(command, 0)
        with patch.dict(os.environ, {"LSB_JOBID": "123", "GU_BATCH_ID": "stale"}), patch.object(hpc.subprocess, "run", side_effect=scientific_stub):
            hpc.execute_task(plan, task)
        self.assertEqual(environments[0]["GU_BATCH_ID"], "batch00002")

    def test_retry_indices(self):
        self.assertEqual(hpc.parse_indices("1,3-5", 9), [1, 3, 4, 5])
        self.assertEqual(hpc.parse_indices("none", 9), [])
        for value in ["0", "10", "2-1", "x"]:
            with self.assertRaises(ValueError):
                hpc.parse_indices(value, 9)

    def test_worker_success_requires_real_canonical_result(self):
        plan = self.make_plan()
        task = plan["tasks"][0]
        hpc.atomic_json(self.run / "state" / "prepare.json", {"plan_id": plan["plan_id"], "state": "complete"})
        with patch.dict(os.environ, {"LSB_JOBID": "123", "LSB_DJOB_NUMPROC": "5"}):
            # Stub program exits 0 but writes nothing: must be FAILED, not zero.
            with self.assertRaisesRegex(ValueError, "no canonical"):
                hpc.execute_task(plan, task)
        receipt = hpc.read_json(hpc.task_receipt(plan, task))
        self.assertEqual(receipt["state"], "failed")
        self.assertNotIn("artifacts", receipt)

    def test_worker_nonzero_exit_propagated(self):
        (self.root / "gu.sh").write_text("#!/bin/bash\nexit 7\n")
        plan = self.make_plan()
        task = plan["tasks"][0]
        hpc.atomic_json(self.run / "state" / "prepare.json", {"plan_id": plan["plan_id"], "state": "complete"})
        with patch.dict(os.environ, {"LSB_JOBID": "123", "LSB_DJOB_NUMPROC": "5"}):
            with self.assertRaises(subprocess.CalledProcessError) as failure:
                hpc.execute_task(plan, task)
        self.assertEqual(failure.exception.returncode, 7)
        self.assertEqual(hpc.read_json(hpc.task_receipt(plan, task))["state"], "failed")

    def test_stale_prepare_receipt_not_reused(self):
        plan = self.make_plan()
        receipt = self.run / "state" / "prepare.json"
        hpc.atomic_json(receipt, {"plan_id": plan["plan_id"], "state": "complete"})
        with patch.dict(os.environ, {"LSB_JOBID": "123", "LSB_DJOB_NUMPROC": "1"}), patch.object(hpc, "prepare", side_effect=ValueError("prepare failed")):
            with self.assertRaisesRegex(ValueError, "prepare failed"):
                hpc.execute(argparse.Namespace(run_dir=str(self.run), stage="prepare"))
        self.assertFalse(receipt.exists())

    def test_assembled_arg_has_verified_source_contract(self):
        plan = self.make_plan("arg")
        for task in plan["tasks"]:
            c = task["chromosomes"][0]
            path = Path(task["output"]) / "arg" / "trace" / "threads"
            path.mkdir(parents=True)
            (path / f"chr{c}.trees").write_bytes(b"synthetic tree bytes for receipt-only test")
            (path / f"chr{c}.sample_map.tsv").write_text("tree_node_id\tsample\thaplotype\n0\ta\t1\n")
            (path / "ARG_TRACE_BUILD.tsv").write_text(f"method\tthreads\nformat\ttrace\nbuild\tGRCh37\nchromosomes\t{c}\n")
            self.seal(plan, task, hpc.artifact_candidates(task))
        hpc.assemble_arg(plan)
        trace = {"grch": "37", "chromosomes": ["1", "2", "X"]}
        self.assertTrue(hpc.validate_arg_assembly(trace, self.run / "arg-assembled"))
        self.assertFalse(hpc.validate_arg_assembly(trace, self.base / "ordinary-native-arg"))
        (self.run / "arg-assembled" / "chr1.sample_map.tsv").write_text("changed sample map")
        with self.assertRaisesRegex(ValueError, "outputs_changed"):
            hpc.validate_arg_assembly(trace, self.run / "arg-assembled")


if __name__ == "__main__":
    unittest.main()
