"""Offline integration checks against an installed, real Kraken 2.

Run with METAG_KRAKEN2=/path/to/kraken2 python3 -m unittest discover
    -s f/global_tests -p test_global_native.py -v

The fixtures are synthetic and small. They exercise the native interface and
publication rules; they are not a throughput or biological-accuracy benchmark.
No reference download, GPU, or production dataset is used.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

FOLDER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FOLDER))

import global_query
from library import file_state


def _run(command, *, input=None):
    result = subprocess.run(list(map(str, command)), input=input,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=120)
    if result.returncode:
        raise AssertionError(f"Native command exited {result.returncode}: {command}\n"
                             + result.stderr.decode(errors="replace"))
    return result


def _taxonomy(folder):
    folder.mkdir(parents=True, exist_ok=True)
    rows = [
        (1, 1, "no rank", "root"),
        (2, 1, "superkingdom", "Bacteria"),
        (100, 2, "genus", "Synthetic genus"),
        (101, 100, "species", "Synthetic bacterium A"),
        (102, 101, "strain", "Synthetic bacterium A strain 1"),
        (103, 100, "species", "Synthetic bacterium B"),
        (104, 103, "strain", "Synthetic bacterium B strain 2"),
        (2759, 1, "superkingdom", "Eukaryota"),
        (4751, 2759, "kingdom", "Fungi"),
        (201, 4751, "species", "Synthetic fungus"),
        (202, 201, "strain", "Synthetic fungus strain 1"),
        (10239, 1, "superkingdom", "Viruses"),
        (301, 10239, "species", "Synthetic virus"),
        (302, 301, "isolate", "Synthetic virus isolate 1"),
    ]
    (folder/"nodes.dmp").write_text("".join(
        f"{tid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
        for tid, parent, rank, _ in rows))
    (folder/"names.dmp").write_text("".join(
        f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n"
        for tid, _, _, name in rows))
    (folder/"merged.dmp").write_text("90002\t|\t102\t|\n")


def _sequences():
    rng = random.Random(20261030)
    dna = lambda n: "".join(rng.choices("ACGT", k=n))
    shared = dna(512)
    unique = [dna(4096) for _ in range(4)]
    refs = [(102, shared+"N"*35+unique[0]),
            (104, shared+"N"*35+unique[1]),
            (202, unique[2]), (302, unique[3])]
    # Shared sequence calls the genus; all N calls unclassified.
    templates = unique + [shared, "N"*512]
    return refs, templates


class NativeGlobalIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.environ.get("METAG_KRAKEN2") or shutil.which("kraken2")
        if not path:
            raise unittest.SkipTest("Set METAG_KRAKEN2 to an installed Kraken 2 wrapper")
        cls.executable = Path(path).resolve(strict=True)
        cls.native = cls.executable.parent/"classify"
        cls.builder = cls.executable.parent/"kraken2-build"
        if not all(p.is_file() and os.access(p, os.X_OK)
                   for p in (cls.executable, cls.native, cls.builder)):
            raise AssertionError("METAG_KRAKEN2 must name a complete executable installation")
        cls.temp = tempfile.TemporaryDirectory(prefix="metag-global-native-")
        cls.addClassCleanup(cls.temp.cleanup)
        cls.root = Path(cls.temp.name)
        cls.refs, cls.templates = _sequences()
        cls.reference = cls.root/"references.fna"
        cls.reference.write_text("".join(
            f">fixture_{i}|kraken:taxid|{tid}\n{seq}\n"
            for i, (tid, seq) in enumerate(cls.refs, 1)))
        cls.full = cls._build("full", 1 << 20)
        cls.sampled = cls._build("sampled", 4096)
        cls.native_version = _run([cls.executable, "--version"]).stdout.decode().strip()
        print("\nReal native fixture:", cls.native_version, flush=True)

    @classmethod
    def _build(cls, name, max_bytes):
        db = cls.root/name
        _taxonomy(db/"taxonomy")
        _run([cls.builder, "--add-to-library", cls.reference,
              "--no-masking", "--db", db])
        built = _run([cls.builder, "--build", "--db", db,
                      "--threads", "1", "--max-db-size", str(max_bytes)])
        (cls.root/(name+".build.log")).write_bytes(built.stdout+built.stderr)
        cls_states = {str(db/name): file_state(db/name)
                      for name in ("hash.k2d", "taxo.k2d", "opts.k2d")}
        return {"path": str(db), "taxonomy_path": str(db/"taxonomy"),
                "outputs": cls_states}

    def setUp(self):
        self.work = self.root/self._testMethodName
        self.work.mkdir()
        self.args = SimpleNamespace(threads=2, confidence=0.1,
                                    minimum_hit_groups=2, minimum_base_quality=0,
                                    min_free_gib=0, memory_reserve_gib=0,
                                    keep_work=True, work_root=self.work/"work",
                                    benchmark_fragments=17)

    def _reads(self, count, *, paired=True):
        left, right, expected = [], [], []
        normalized = (101, 103, 201, 301, 100, 0)
        for index in range(count):
            category = index % 6
            template = self.templates[category]
            offset = (index*13) % (len(template)-240)
            nl = b"\r\n" if index % 7 == 0 else b"\n"
            quality = bytes(33+(index+j)%41 for j in range(135))
            records = []
            for mate, start, length in ((1, offset, 120), (2, offset+83, 135)):
                # Duplicate IDs are intentional; ordinal position is identity.
                header = f"@duplicate/{mate} original-index={index}".encode()+nl
                seq = template[start:start+length].encode()+nl
                plus = f"+duplicate/{mate} original-index={index}".encode()+nl
                records.append((header, seq, plus, quality[:length]+nl))
            left.append(b"".join(records[0]))
            right.append(b"".join(records[1]))
            expected.append(normalized[category])
        r1 = self.work/"input_R1.fastq"
        r2 = self.work/"input_R2.fastq.gz"
        r1.write_bytes(b"".join(left))
        with gzip.open(r2, "wb") as handle:
            handle.write(b"".join(right))
        sample = {"name": "native_fixture", "r1": r1,
                  "r2": r2 if paired else None, "target": self.work/"split"}
        return sample, expected, left, right

    @staticmethod
    def _taxids(query):
        return [row[0] for row in struct.iter_unpack("<Q", Path(query["taxids"]).read_bytes())]

    def test_native_numeric_ids_cross_thread_blocks_and_species_normalization(self):
        template = self.templates[0]
        payload = (f"@123\n{template[:120]}\n+\n{'I'*120}\n"
                   f"@123\n{template[83:218]}\n+\n{'I'*135}\n").encode()
        direct = _run(global_query.classifier_command(self.native, self.full["path"],
                                                      self.args, True), input=payload)
        fields = direct.stdout.rstrip(b"\n").split(b"\t")
        self.assertEqual(fields[:4], [b"C", b"123", b"102", b"120|135"])
        sample, expected, _, _ = self._reads(20_013)
        result = global_query.query_once(sample, self.full, self.executable,
                                         self.args, self.work/"query")
        self.assertEqual(result["fragments"], len(expected))
        self.assertEqual(result["reads"], 2*len(expected))
        self.assertEqual(result["compact_bytes"], 8*len(expected))
        self.assertEqual(self._taxids(result), expected)
        self.assertIn("-S", result["command"])
        self.assertFalse((self.work/"query"/"lengths.partial.bin").exists())
        self.assertFalse((self.work/"query"/"seen.partial.bin").exists())

    def test_quick_sampled_database_is_thread_consistent_and_keeps_all_pairs(self):
        sample, _, _, _ = self._reads(20_013)
        self.args.quick = True
        self.args.confidence = 0
        outputs = []
        for threads in (1, 4):
            self.args.threads = threads
            query = global_query.query_once(sample, self.sampled, self.executable,
                                            self.args, self.work/f"quick-{threads}")
            self.assertEqual(query["fragments"], 20_013)
            self.assertEqual(query["reads"], 40_026)
            self.assertIn("-q", query["command"])
            outputs.append(self._taxids(query))
        self.assertEqual(outputs[0], outputs[1])
        self.assertTrue(set(outputs[0]) <= {0, 100, 101, 103, 201, 301})

    def test_split_counts_original_fastq_bytes_and_completed_reuse(self):
        sample, expected, left, right = self._reads(61)
        output = global_query.run_sample(sample, self.full, self.executable,
                                         self.args)
        marker = json.loads((output/"split.done").read_text())
        self.assertEqual((marker["fragments"], marker["reads"]), (61, 122))
        self.assertEqual(marker["reference_evidence"], "not_recorded")
        with (output/"counts.tsv").open() as handle:
            counts = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(sum(int(row["fragments"]) for row in counts), 61)
        self.assertEqual(sum(int(row["reads"]) for row in counts), 122)
        self.assertEqual({int(row["assigned_taxid"]): int(row["fragments"])
                          for row in counts}, dict(Counter(expected)))
        expected_bins = defaultdict(lambda: [bytearray(), bytearray()])
        for tid, a, b in zip(expected, left, right):
            expected_bins[tid][0].extend(a)
            expected_bins[tid][1].extend(b)
        folders = {101: "bacteria/101_Synthetic_bacterium_A",
                   103: "bacteria/103_Synthetic_bacterium_B",
                   201: "fungi/201_Synthetic_fungus",
                   301: "virus/301_Synthetic_virus",
                   100: "split.bad/bacteria/no_species_in_taxonomic_lineage",
                   0: "split.bad/unclassified/unclassified"}
        for tid, mates in expected_bins.items():
            for mate in (1, 2):
                with gzip.open(output/folders[tid]/f"PE.R{mate}.fastq.gz", "rb") as handle:
                    self.assertEqual(handle.read(), bytes(mates[mate-1]))
        with gzip.open(output/"read-assignments.tsv.gz", "rt") as handle:
            assignments = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(assignments), 122)
        self.assertEqual([int(row["fragment_index"]) for row in assignments],
                         [index for index in range(1, 62) for _ in range(2)])
        with gzip.open(output/"read-reference-hits.tsv.gz", "rt") as handle:
            self.assertEqual(len(list(csv.reader(handle, delimiter="\t"))), 1)
        state = file_state(output/"split.done")
        self.assertEqual(global_query.run_sample(sample, self.full, self.executable,
                                                self.args), output)
        self.assertEqual(file_state(output/"split.done"), state)

    def test_benchmark_prefix_does_not_publish_species_output(self):
        sample, _, _, _ = self._reads(61)
        work = global_query.benchmark(sample, self.full, self.executable, self.args)
        report = json.loads((work/"benchmark.json").read_text())
        self.assertEqual(report["query"]["fragments"], 17)
        self.assertEqual(report["query"]["validated_scope"], "prefix")
        self.assertFalse(report["species_output_created"])
        self.assertTrue(report["includes_database_load"])
        self.assertFalse(sample["target"].exists())
        self.assertFalse(list(work.rglob("*.fastq.gz")))
        self.assertFalse(list(work.rglob("split.done")))

    def test_single_end_native_mode_and_counts(self):
        sample, expected, _, _ = self._reads(31, paired=False)
        result = global_query.query_once(sample, self.full, self.executable,
                                         self.args, self.work/"query")
        self.assertEqual(result["reads"], 31)
        self.assertFalse(result["paired"])
        self.assertNotIn("-S", result["command"])
        self.assertEqual(self._taxids(result), expected)

    def test_sampled_database_native_regression_and_thread_consistency(self):
        # Force hash subsampling even for this tiny reference library. This
        # exercises the skipped-lookup/repeat branch fixed upstream in
        # c6fc2bf02cd6a52d1dd6d472b1b8a59cca5d94e2 (included in 2.17.2).
        sampled_path = Path(self.sampled["path"])
        inspected = _run([self.executable.parent/"kraken2-inspect", "--db", sampled_path])
        diagnostic = (inspected.stdout+inspected.stderr).decode(errors="replace")
        self.assertIn("Min clear hash value", diagnostic)
        line = next(line for line in diagnostic.splitlines()
                    if "Min clear hash value" in line)
        self.assertGreater(int(line.rsplit("=", 1)[-1].strip()), 0)
        sample, expected, _, _ = self._reads(1201)
        self.args.threads = 1
        first = global_query.query_once(sample, self.sampled, self.executable,
                                        self.args, self.work/"one_thread")
        self.args.threads = 2
        second = global_query.query_once(sample, self.sampled, self.executable,
                                         self.args, self.work/"two_threads")
        actual = self._taxids(first)
        self.assertEqual(actual, self._taxids(second))
        self.assertEqual(len(actual), len(expected))
        self.assertTrue(any(tid != 0 for tid in actual))
        self.assertTrue(any(tid == 0 and truth != 0 for tid, truth in zip(actual, expected)))
        for tid, truth in zip(actual, expected):
            self.assertIn(tid, (0, truth))

    def test_global_builder_real_native_add_mask_build_inspect_and_query(self):
        """Exercise the public builder against the same actual native programs."""
        import global_build
        import library

        lib, raw = self.work/"library", self.work/"raw"
        inventory = lib/"reference/inventory.tsv"
        inventory.parent.mkdir(parents=True)
        (raw/"microbes").mkdir(parents=True)
        rows, original, states = [], {}, {}
        for index, (taxid, seq) in enumerate(self.refs, 1):
            suffix = ".fna.gz" if index == 2 else ".fna"
            relative = f"microbes/reference{index}{suffix}"
            filename = raw/relative
            content = f">original_contig_{index} retained source header\n{seq}\n".encode()
            if suffix.endswith(".gz"):
                with gzip.open(filename, "wb") as handle:
                    handle.write(content)
            else:
                filename.write_bytes(content)
            original[filename] = filename.read_bytes()
            info = filename.stat()
            states[relative] = [info.st_size, info.st_mtime_ns, info.st_ctime_ns,
                                info.st_dev, info.st_ino]
            rows.append({"relative_path": relative, "filename": filename.name,
                         "taxid": str(taxid)})
        with inventory.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("relative_path", "filename", "taxid"),
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        info = inventory.stat()
        (inventory.parent/"check-ref.receipt.json").write_text(json.dumps({
            "status": "complete", "references": states,
            "outputs": {str(inventory): [info.st_mode, info.st_size, info.st_mtime_ns,
                                         info.st_ctime_ns, info.st_dev, info.st_ino]}}))
        _taxonomy(lib/"metacache/taxonomy")
        args = SimpleNamespace(db=self.work/"database", db_gib=(1 << 20)/2**30,
                               min_free_gib=0, memory_reserve_gib=0, threads=1,
                               chunk_bases=4096, no_masking=False)
        # Only fixture paths and PATH are configured here. Native add, masker,
        # builder, inspector, memory/disk checks and classification are real.
        with ExitStack() as stack:
            for module in (global_build, library):
                stack.enter_context(patch.object(module, "DIR_LIB", lib))
                stack.enter_context(patch.object(module, "DIR_REF", raw))
            stack.enter_context(patch.dict(os.environ, {
                "PATH": str(self.executable.parent)+os.pathsep+os.environ.get("PATH", "")}))
            receipt = global_build.build_database(args)
            self.assertEqual(receipt["status"], "complete")
            self.assertEqual(receipt["reference"]["active_references"], 4)
            self.assertEqual(receipt["library"]["sequence_records"], 4)
            self.assertEqual(receipt["library"]["original_bases"], sum(len(s) for _, s in self.refs))
            self.assertGreater(receipt["library"]["overlap_bases"], 0)
            self.assertTrue(receipt["native_inspect_passed"])
            self.assertTrue(receipt["parameters"]["masking"])
            snapshot = global_build.database_snapshot(args.db)
            self.assertEqual(snapshot, receipt)
            native_before = {name: file_state(args.db/name)
                             for name in ("hash.k2d", "taxo.k2d", "opts.k2d")}
            self.assertEqual(global_build.build_database(args), receipt)
            self.assertEqual({name: file_state(args.db/name) for name in native_before}, native_before)
            sample, expected, _, _ = self._reads(61)
            result = global_query.query_once(sample, snapshot, self.executable,
                                             self.args, self.work/"built_query")
            self.assertEqual(self._taxids(result), expected)
            for filename, content in original.items():
                self.assertEqual(filename.read_bytes(), content)


if __name__ == "__main__":
    unittest.main()
