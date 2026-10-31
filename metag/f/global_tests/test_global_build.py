"""Focused input, chunk, and interrupted-publication tests for global build.

Native commands are simulated here to force failures at specific boundaries.
test_global_native.py separately exercises the actual installed Kraken2 tools.
"""
from contextlib import ExitStack
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import global_build as build
import library


def fasta_records(payload):
    records = []
    for line in payload.splitlines():
        if line.startswith(b">"):
            records.append([line[1:], bytearray()])
        else:
            records[-1][1].extend(line)
    return [(header, bytes(sequence)) for header, sequence in records]


class BuildFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="metag-global-build-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.lib, self.ref, self.db = self.root/"lib", self.root/"refs", self.root/"db"
        (self.ref/"microbes").mkdir(parents=True)
        inventory = self.lib/"reference/inventory.tsv"
        inventory.parent.mkdir(parents=True)
        rng = random.Random(129)
        self.sequences = ["".join(rng.choices("ACGT", k=n)).encode() for n in (40, 60, 180, 25)]
        rows, self.original, states = [], {}, {}
        for index, sequences in enumerate((self.sequences[:2], self.sequences[2:]), 1):
            relative = f"microbes/reference{index}.fna" + (".gz" if index == 2 else "")
            path = self.ref/relative
            payload = b"".join(b">duplicate source name\n"+seq+b"\n" for seq in sequences)
            if index == 2:
                with gzip.open(path, "wb") as handle:
                    handle.write(payload)
            else:
                path.write_bytes(payload)
            self.original[path] = path.read_bytes()
            info = path.stat()
            states[relative] = [info.st_size, info.st_mtime_ns, info.st_ctime_ns, info.st_dev, info.st_ino]
            rows.append({"relative_path": relative, "filename": path.name, "taxid": str(index+1)})
        self.excluded = self.ref/"microbes/excluded.fna"
        self.excluded.write_text(">must_not_be_included\n"+"Z"*200+"\n")
        with inventory.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("relative_path", "filename", "taxid"),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        info = inventory.stat()
        self.validation = {"status": "complete", "references": states,
                           "outputs": {str(inventory): [info.st_mode, info.st_size, info.st_mtime_ns,
                                                        info.st_ctime_ns, info.st_dev, info.st_ino]}}
        self.validation_path = inventory.parent/"check-ref.receipt.json"
        self.validation_path.write_text(json.dumps(self.validation))
        taxonomy = self.lib/"metacache/taxonomy"
        taxonomy.mkdir(parents=True)
        (taxonomy/"nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tspecies\t|\n3\t|\t1\t|\tspecies\t|\n")
        (taxonomy/"names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        (taxonomy/"merged.dmp").write_text("99\t|\t2\t|\n")
        self.stack = self.enterContext(ExitStack())
        for module in (build, library):
            self.stack.enter_context(patch.object(module, "DIR_LIB", self.lib))
            self.stack.enter_context(patch.object(module, "DIR_REF", self.ref))
        self.stack.enter_context(patch.object(build, "DIR_RAW", self.root/"reads"))
        self.stack.enter_context(patch.object(build, "ROOT", self.root/"code"))
        self.stack.enter_context(patch.object(build, "space_guard"))
        self.memory = self.stack.enter_context(patch.object(build, "require_memory", return_value={}))
        self.args = SimpleNamespace(db=self.db, db_gib=1/1024, min_free_gib=0,
                                    memory_reserve_gib=0, threads=2, chunk_bases=100,
                                    no_masking=False)
        self.programs = {"kraken2-build": "/fixture/kraken2-build",
                         "kraken2-inspect": "/fixture/kraken2-inspect",
                         "version": "Kraken version 2.17.1", "wrapper_sha256": "wrapper",
                         "native_sha256": {"build_db": "native"}}
        self.stack.enter_context(patch.object(build, "_software", return_value=self.programs))
        self.commands = []
        self.fail_add_number = None
        self.fail_inspect_once = False
        self.adds = 0
        self.native = self.stack.enter_context(patch.object(build, "_run", side_effect=self.fake_native))

    def fake_native(self, command, log, db, reserve):
        self.commands.append(list(command))
        folder = Path(command[command.index("--db")+1])
        if "--add-to-library" in command:
            self.adds += 1
            source = Path(command[command.index("--add-to-library")+1])
            added = folder/"library/added"
            added.mkdir(parents=True)
            shutil.copyfile(source, added/"native_random.fna")
            records = fasta_records(source.read_bytes())
            (added/"prelim_map_random.txt").write_text("".join(
                f"TAXID\t{header.decode()}\t{header.rsplit(b'|', 1)[-1].decode()}\n"
                for header, _ in records))
            if self.adds == self.fail_add_number:
                raise RuntimeError("simulated interrupted native add")
        elif "--build" in command:
            for name in build.NATIVE_FILES:
                (folder/name).write_bytes(("fake native "+name).encode())
        elif command[0].endswith("kraken2-inspect"):
            if self.fail_inspect_once:
                self.fail_inspect_once = False
                raise RuntimeError("simulated invalid native index")
        else:
            raise AssertionError(command)
        Path(log).write_text("simulated native command\n")
        return 0.01

    def test_only_validated_active_references_and_taxids_are_used(self):
        rows, source = build._verified_inputs()
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["taxid"] for row in rows}, {2, 3})
        self.assertNotIn(str(self.excluded), [row["path"] for row in rows])
        self.assertEqual(source["active_references"], 2)
        self.assertTrue(source["excluded_references_are_not_added"])
        changed = next(iter(self.original))
        changed.write_bytes(changed.read_bytes()+b"A")
        with self.assertRaisesRegex(RuntimeError, "Reference differs"):
            build._verified_inputs()

    def test_changed_inventory_is_rejected(self):
        inventory = self.lib/"reference/inventory.tsv"
        inventory.write_text(inventory.read_text().replace("\t2\n", "\t99\n"))
        with self.assertRaisesRegex(RuntimeError, "inventory differs"):
            build._verified_inputs()

    def test_missing_taxid_is_rejected_after_valid_inventory_receipt(self):
        taxonomy = self.lib/"metacache/taxonomy/nodes.dmp"
        taxonomy.write_text("1\t|\t1\t|\tno rank\t|\n")
        with self.assertRaisesRegex(RuntimeError, "taxids absent"):
            build._verified_inputs()

    def test_chunks_preserve_every_35mer_and_do_not_split_short_records(self):
        from collections import Counter
        rows, _ = build._verified_inputs()
        work = self.root/"work"
        work.mkdir()
        captured, stats = [], []
        for chunk in build._formatted_chunks(rows, work, 100, 0):
            stats.append(dict(chunk))
            captured.extend(fasta_records(chunk["path"].read_bytes()))
            self.assertLessEqual(chunk["bases"], 100)
        self.assertEqual(sum(item["original_bases"] for item in stats), sum(map(len, self.sequences)))
        self.assertEqual(sum(item["sequence_records"] for item in stats), 4)
        self.assertEqual(len({header for header, _ in captured}), len(captured))
        self.assertEqual(sum(len(seq) for _, seq in captured)-sum(map(len, self.sequences)), 2*34)
        self.assertEqual(sorted(len(seq) for _, seq in captured if len(seq) < 100), [25, 40, 48, 60])
        expected = Counter(seq[i:i+35] for seq in self.sequences for i in range(len(seq)-34))
        actual = Counter(seq[i:i+35] for _, seq in captured for i in range(len(seq)-34))
        self.assertEqual(actual, expected)
        self.assertFalse((work/"current-chunk.fna").exists())

    def test_completed_database_reuses_native_work_and_keeps_source_unchanged(self):
        receipt = build.build_database(self.args)
        calls = len(self.commands)
        self.assertEqual(build.build_database(self.args), receipt)
        self.assertEqual(len(self.commands), calls)
        self.assertEqual(build.database_snapshot(self.db), receipt)
        self.assertTrue(receipt["parameters"]["masking"])
        self.assertEqual(receipt["library"]["original_bases"], 305)
        retained = b"".join(path.read_bytes() for path in (self.db/"library").rglob("*.fna"))
        self.assertNotIn(b"Z", retained)
        self.assertFalse((self.db/".global-work/current-chunk.fna").exists())
        for path, payload in self.original.items():
            self.assertEqual(path.read_bytes(), payload)
        adds = [cmd for cmd in self.commands if "--add-to-library" in cmd]
        self.assertTrue(adds)
        self.assertTrue(all("--no-masking" not in cmd for cmd in adds))
        self.assertTrue(all(str(self.db/"taxonomy"/name) in receipt["outputs"] for name in build.TAXONOMY_FILES))

    def test_no_masking_is_only_an_explicit_option(self):
        self.args.no_masking = True
        receipt = build.build_database(self.args)
        self.assertFalse(receipt["parameters"]["masking"])
        self.assertTrue(all("--no-masking" in cmd for cmd in self.commands if "--add-to-library" in cmd))
        with patch.dict(os.environ, {"KRAKEN2_MASK_LC": "0", "KRAKEN2_KMER_LEN": "17"}):
            env = build._environment()
        self.assertEqual(env["KRAKEN2_MASK_LC"], "1")
        self.assertNotIn("KRAKEN2_KMER_LEN", env)

    def test_interrupted_add_reuses_committed_chunk_without_duplicate_append(self):
        self.fail_add_number = 2
        with self.assertRaisesRegex(RuntimeError, "interrupted native add"):
            build.build_database(self.args)
        first = self.db/"library/chunk00000001/sequences.fna"
        first_state = library.file_state(first)
        self.assertFalse((self.db/"build.json").exists())
        self.assertTrue((self.db/".global-work/adding").exists())
        receipt = build.build_database(self.args)
        self.assertEqual(library.file_state(first), first_state)
        self.assertEqual(self.adds, receipt["library"]["chunks"]+1)
        self.assertEqual(receipt["library"]["original_bases"], 305)
        self.assertFalse((self.db/".global-work/adding").exists())

    def test_changed_committed_chunk_is_never_treated_as_complete(self):
        self.fail_add_number = 2
        with self.assertRaises(RuntimeError):
            build.build_database(self.args)
        first = self.db/"library/chunk00000001/sequences.fna"
        first.write_bytes(b"corrupt")
        with self.assertRaisesRegex(RuntimeError, "incomplete/changed"):
            build.build_database(self.args)
        self.assertFalse((self.db/"build.json").exists())

    def test_inspect_failure_keeps_chunks_but_rebuilds_native_files_on_retry(self):
        self.fail_inspect_once = True
        with self.assertRaisesRegex(RuntimeError, "invalid native index"):
            build.build_database(self.args)
        adds = self.adds
        self.assertFalse((self.db/"build.json").exists())
        self.assertTrue((self.db/".global-work/library-ready.json").exists())
        self.assertTrue(all(not (self.db/name).exists() for name in build.NATIVE_FILES))
        receipt = build.build_database(self.args)
        self.assertEqual(self.adds, adds)
        self.assertEqual(receipt["status"], "complete")
        self.assertEqual(sum("--build" in cmd for cmd in self.commands), 2)

    def test_hard_kill_with_all_native_files_rebuilds_uncommitted_hash(self):
        self.fail_inspect_once = True
        with self.assertRaises(RuntimeError):
            build.build_database(self.args)
        adds = self.adds
        for name in build.NATIVE_FILES:
            (self.db/name).write_bytes(b"present but not committed after OOM")
            (self.db/(name + ".tmp")).write_bytes(b"incomplete native scratch")
        def native(command, log, db, reserve):
            if "--build" in command:
                for name in build.NATIVE_FILES:
                    self.assertFalse((self.db/name).exists())
                    self.assertFalse((self.db/(name + ".tmp")).exists())
            return self.fake_native(command, log, db, reserve)
        self.native.side_effect = native
        result = build.build_database(self.args)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(self.adds, adds)

    def test_snapshot_detects_changed_frozen_taxonomy(self):
        build.build_database(self.args)
        taxonomy = self.db/"taxonomy/nodes.dmp"
        taxonomy.write_bytes(taxonomy.read_bytes()+b"4\t|\t1\t|\n")
        with self.assertRaisesRegex(RuntimeError, "No unchanged, completed"):
            build.database_snapshot(self.db)

    def test_memory_failure_happens_before_output_mutation(self):
        self.memory.side_effect = RuntimeError("insufficient effective memory")
        with self.assertRaisesRegex(RuntimeError, "insufficient effective memory"):
            build.build_database(self.args)
        self.assertFalse(self.db.exists())
        self.assertFalse(self.commands)

    def test_plan_reports_insufficient_resources_without_native_execution(self):
        with patch.object(build, "memory_status", return_value={"effective_available_bytes": 1}), \
                patch.object(build, "_disk_available", return_value=1):
            result = build.plan_build(db=self.db)
        self.assertFalse(result["memory_fits"])
        self.assertFalse(result["initial_disk_fits"])
        self.assertEqual(result["parameters"]["hash_cap_bytes"], 64*2**30)
        self.assertFalse(self.db.exists())
        self.assertFalse(self.commands)

    def test_unowned_directory_is_not_adopted(self):
        self.db.mkdir()
        foreign = self.db/"unrelated.txt"
        foreign.write_text("retain this")
        with self.assertRaisesRegex(RuntimeError, "nonempty, unowned"):
            build.build_database(self.args)
        self.assertEqual(foreign.read_text(), "retain this")
        self.assertFalse(self.commands)


if __name__ == "__main__":
    unittest.main()
