"""Public CLI selection, provenance and resource admission boundaries."""
from contextlib import redirect_stdout, redirect_stderr
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import global_pipeline as cli
import global_resources as resources


def install_receipt(folder, reported="Kraken version 2.17.1"):
    receipt = {"source_tag": cli.KRAKEN_VERSION, "source_commit": cli.KRAKEN_COMMIT,
               "reported_version": reported,
               "files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in folder.iterdir() if p.is_file()}}
    path = folder.parent / "install.json"
    path.write_text(json.dumps(receipt))
    return path


class CliTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="metag-global-cli-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.folder = self.root / "installation/bin"
        self.folder.mkdir(parents=True)
        for name in cli.KRAKEN_REQUIRED_FILES:
            path = self.folder / name
            path.write_text("#!/bin/sh\nexit 0\n")
            path.chmod(0o755)
        (self.folder / "kraken2").write_text("#!/bin/sh\nprintf 'Kraken version 2.17.1\\n'\n")
        self.receipt_path = install_receipt(self.folder)
        self.args = SimpleNamespace(kraken_dir=self.folder)
        self.enterContext(patch.object(cli, "ensure_memory_scope", return_value=None))
        self.enterContext(patch.dict(os.environ, {
            "METAG_TOOLS": str(self.root / "tools"), "METAG_KRAKEN2_DIR": ""}))

    def test_pinned_tag_accepts_official_old_version_label_and_records_both(self):
        result = cli.configure_tools(self.args, required=True)
        self.assertEqual(result, self.folder / "kraken2")
        metadata = self.args.kraken_installation
        self.assertEqual(metadata["source_tag"], "2.17.2")
        self.assertEqual(metadata["reported_version"], "2.17.1")
        self.assertTrue(metadata["version_label_mismatch"])
        self.assertEqual(metadata["source_commit"], cli.KRAKEN_COMMIT)
        self.assertEqual(metadata["install_receipt_sha256"],
                         hashlib.sha256(self.receipt_path.read_bytes()).hexdigest())

    def test_plain_2171_without_pinned_receipt_is_rejected(self):
        self.receipt_path.unlink()
        with patch.object(cli.subprocess, "run") as probe:
            with self.assertRaisesRegex(RuntimeError, "receipt"):
                cli.configure_tools(self.args, required=True)
        probe.assert_not_called()

    def test_modified_binary_is_rejected_before_wrapper_execution(self):
        (self.folder / "classify").write_text("modified native binary\n")
        with patch.object(cli.subprocess, "run") as probe:
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                cli.configure_tools(self.args, required=True)
        probe.assert_not_called()

    def test_wrong_source_missing_manifest_entry_and_extra_file_are_rejected(self):
        original = json.loads(self.receipt_path.read_text())
        for mode in ("commit", "manifest", "extra"):
            with self.subTest(mode=mode):
                receipt = json.loads(json.dumps(original))
                if mode == "commit":
                    receipt["source_commit"] = "0" * 40
                elif mode == "manifest":
                    del receipt["files"]["build_db"]
                else:
                    (self.folder / "unrecorded-helper").write_text("unrecorded\n")
                self.receipt_path.write_text(json.dumps(receipt))
                with self.assertRaises(RuntimeError):
                    cli.configure_tools(self.args, required=True)

    def test_version_receipt_and_executable_permissions_must_match(self):
        install_receipt(self.folder, reported="2.17.2")
        with self.assertRaisesRegex(RuntimeError, "receipt reports"):
            cli.configure_tools(self.args, required=True)
        install_receipt(self.folder)
        (self.folder / "build_db").chmod(0o644)
        with self.assertRaisesRegex(RuntimeError, "Incomplete"):
            cli.configure_tools(self.args, required=True)

    def test_explicit_incomplete_directory_cannot_fall_back_to_path(self):
        other = self.root / "empty-bin"
        other.mkdir()
        with patch.dict(os.environ, {"PATH": str(self.folder) + os.pathsep + os.environ["PATH"]}):
            with self.assertRaisesRegex(RuntimeError, "missing in"):
                cli.configure_tools(SimpleNamespace(kraken_dir=other), required=True)

    def test_tool_plan_does_not_execute_wrapper_or_change_path(self):
        old_path = os.environ["PATH"]
        with patch.object(cli.subprocess, "run", side_effect=AssertionError("plan executed a tool")):
            report = cli.configure_tools(self.args, required=False)
        self.assertEqual(report["executable"], str(self.folder / "kraken2"))
        self.assertEqual(os.environ["PATH"], old_path)
        self.assertEqual(report["required_source_commit"], cli.KRAKEN_COMMIT)

    def test_cli_defaults_and_invalid_resource_or_classifier_values(self):
        args = cli.parser().parse_args(["split"])
        self.assertEqual(args.output_suffix, "kraken2")
        self.assertEqual(args.work_root, cli.DIR_OUT / "runs_global")
        self.assertEqual((args.confidence, args.minimum_hit_groups, args.minimum_base_quality), (0, 2, 0))
        self.assertEqual(args.memory_reserve_gib, 8)
        self.assertFalse(args.keep_work)
        self.assertEqual(cli.parser().parse_args(["benchmark"]).benchmark_fragments, 1_000_000)
        for option, value in (("--confidence", "nan"), ("--confidence", "1.01"),
                              ("--min-free-gib", "-1"), ("--threads", "0"),
                              ("--minimum-hit-groups", "0"), ("--minimum-base-quality", "94")):
            with self.subTest(option=option, value=value), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    cli.parser().parse_args(["split", option, value])
                self.assertEqual(error.exception.code, 2)

    def test_split_dispatch_passes_matching_installation_and_separate_targets(self):
        import global_build
        import global_query
        sample = {"name": "Ready__run__kraken2", "target": self.root / "split/Ready/run__kraken2"}
        snapshot = {"path": str(self.root / "database"), "identity": "fixture"}
        with patch.object(cli, "assert_reference_policy"), \
                patch.object(cli, "discover_samples", return_value=[sample]), \
                patch.object(global_build, "database_snapshot", return_value=snapshot), \
                patch.object(global_query, "run_sample") as run_sample, redirect_stdout(io.StringIO()):
            status = cli.main(["split", "--kraken-dir", str(self.folder), "--sample", "Ready/run"])
        self.assertEqual(status, 0)
        run_sample.assert_called_once()
        actual_sample, actual_snapshot, executable, args = run_sample.call_args.args
        self.assertIs(actual_sample, sample)
        self.assertIs(actual_snapshot, snapshot)
        self.assertEqual(executable, self.folder / "kraken2")
        self.assertEqual(args.kraken_installation["source_commit"], cli.KRAKEN_COMMIT)
        self.assertEqual(args.work_root.name, "runs_global")

    def test_benchmark_rejects_multiple_samples_before_query(self):
        import global_build
        import global_query
        with patch.object(cli, "assert_reference_policy"), \
                patch.object(cli, "discover_samples", return_value=[{}, {}]), \
                patch.object(global_build, "database_snapshot", return_value={}), \
                patch.object(global_query, "benchmark") as benchmark:
            with self.assertRaisesRegex(RuntimeError, "exactly one sample"):
                cli.main(["benchmark", "--kraken-dir", str(self.folder)])
        benchmark.assert_not_called()

    def test_public_shell_plan_is_read_only_and_uses_requested_chunk_size(self):
        ref, lib, raw = self.root / "refs", self.root / "lib", self.root / "raw"
        reference = ref / "microbes/ref.fna"
        reference.parent.mkdir(parents=True)
        reference.write_text(">reference\n" + "ACGT" * 20 + "\n")
        inventory = lib / "reference/inventory.tsv"
        inventory.parent.mkdir(parents=True)
        inventory.write_text("relative_path\tfilename\ttaxid\nmicrobes/ref.fna\tref.fna\t2\n")
        ref_stat, inv_stat = reference.stat(), inventory.stat()
        (inventory.parent / "check-ref.receipt.json").write_text(json.dumps({
            "status": "complete",
            "references": {"microbes/ref.fna": [ref_stat.st_size, ref_stat.st_mtime_ns,
                          ref_stat.st_ctime_ns, ref_stat.st_dev, ref_stat.st_ino]},
            "outputs": {str(inventory): [inv_stat.st_mode, inv_stat.st_size, inv_stat.st_mtime_ns,
                         inv_stat.st_ctime_ns, inv_stat.st_dev, inv_stat.st_ino]}}))
        taxonomy = lib / "metacache/taxonomy"
        taxonomy.mkdir(parents=True)
        (taxonomy / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tspecies\t|\n")
        (taxonomy / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        (taxonomy / "merged.dmp").write_text("99\t|\t2\t|\n")
        reads = raw / "Ready/run.fastq"
        reads.parent.mkdir(parents=True)
        reads.write_text("@one\nACGT\n+\n!!!!\n")
        old_work = self.root / "output/runs/existing/query.partial"
        old_work.parent.mkdir(parents=True)
        old_work.write_bytes(b"old MetaCache work must remain untouched\n")
        config = self.root / "metag.conf"
        config.write_text("\n".join(f"{name}={path}" for name, path in {
            "DIR_REF": ref, "DIR_LIB": lib, "DIR_RAW": raw,
            "DIR_OUT": self.root / "output", "DIR_SPLIT": self.root / "split"}.items()) +
            "\nKEEP_REFERENCE_FASTA=TRUE\n")
        # A tool probe would create an unexpected file and fail the tree check.
        (self.folder / "kraken2").write_text("#!/bin/sh\ntouch '" + str(self.root / "unexpected-probe") + "'\nexit 99\n")
        def tree():
            return {str(p.relative_to(self.root)): (p.read_bytes(), p.stat().st_mtime_ns)
                    for p in self.root.rglob("*") if p.is_file()}
        before = tree()
        process = subprocess.run(["bash", str(Path(cli.__file__).parents[1] / "3.global.sh"),
                                  "plan", "--kraken-dir", str(self.folder), "--chunk-mib", "3",
                                  "--sample", "Ready/run"], capture_output=True, text=True, timeout=30,
                                 env={**os.environ, "METAG_CONFIG": str(config)})
        self.assertEqual(process.returncode, 0, process.stderr + process.stdout)
        report = json.loads(process.stdout)
        self.assertTrue(report["plan_complete"])
        self.assertEqual(report["samples"]["count"], 1)
        self.assertTrue(report["samples"]["first_sample"]["target"].endswith("run__kraken2"))
        self.assertEqual(report["build"]["parameters"]["chunk_bases"], 3 << 20)
        self.assertEqual(report["build"]["initial_disk_headroom_bytes"],
                         (24 + 20 + 2) * resources.GIB + 4 * (3 << 20))
        self.assertEqual(tree(), before)
        self.assertFalse((lib / "kraken2_global").exists())
        self.assertFalse((self.root / "output/runs_global").exists())
        self.assertFalse((self.root / "split").exists())

    @unittest.skipUnless(os.environ.get("METAG_KRAKEN2"), "set METAG_KRAKEN2 for actual installed native tools")
    def test_actual_official_wrapper_and_complete_installation_receipt(self):
        source = Path(os.environ["METAG_KRAKEN2"]).resolve().parent
        destination = self.root / "actual/bin"
        shutil.copytree(source, destination)
        reported = subprocess.check_output([str(destination / "kraken2"), "--version"], text=True).splitlines()[0]
        install_receipt(destination, reported=reported)
        args = SimpleNamespace(kraken_dir=destination)
        self.assertEqual(cli.configure_tools(args, required=True), destination / "kraken2")
        self.assertEqual(args.kraken_installation["reported_version"], "2.17.1")


class ResourceBoundaryTests(unittest.TestCase):
    def test_exact_available_memory_passes_and_one_byte_less_fails(self):
        required, reserve = 123, resources.GIB
        with patch.object(resources, "memory_status", return_value={"effective_available_bytes": required + reserve}):
            state = resources.require_memory(required, reserve_gib=1)
        self.assertEqual(state["required_bytes"], required)
        self.assertEqual(state["reserve_bytes"], reserve)
        with patch.object(resources, "memory_status", return_value={"effective_available_bytes": required + reserve - 1}):
            with self.assertRaisesRegex(RuntimeError, "Insufficient"):
                resources.require_memory(required, reserve_gib=1)

    def test_v1_unlimited_parent_does_not_hide_finite_child(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            proc, group = root / "proc", root / "cgroup"
            (proc / "self").mkdir(parents=True)
            (proc / "meminfo").write_text("MemTotal: 2097152 kB\nMemAvailable: 1048576 kB\n")
            (proc / "self/cgroup").write_text("5:memory:/job\n")
            child = group / "memory/job"
            child.mkdir(parents=True)
            (child.parent / "memory.limit_in_bytes").write_text(str((1 << 63) - 4096))
            (child / "memory.limit_in_bytes").write_text(str(512 << 20))
            (child / "memory.usage_in_bytes").write_text(str(128 << 20))
            status = resources.memory_status(proc_root=proc, cgroup_root=group)
            self.assertEqual(status["effective_available_bytes"], 384 << 20)
            (child / "memory.usage_in_bytes").unlink()
            with self.assertRaisesRegex(RuntimeError, "finite cgroup"):
                resources.memory_status(proc_root=proc, cgroup_root=group)

    def test_malformed_cgroup_limit_cannot_be_treated_as_unlimited(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            proc, group = root / "proc", root / "cgroup"
            proc.mkdir(); group.mkdir()
            (proc / "meminfo").write_text("MemTotal: 2097152 kB\nMemAvailable: 1048576 kB\n")
            (group / "memory.max").write_text("unreadable-corrupt-limit")
            with self.assertRaisesRegex(RuntimeError, "Invalid cgroup memory value"):
                resources.memory_status(proc_root=proc, cgroup_root=group)


if __name__ == "__main__":
    unittest.main()
