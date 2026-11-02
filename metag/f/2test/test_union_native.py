"""Real CPU classification over two compatible independently built indexes.

Run explicitly (compilation takes time):
  METAG_RUN_NATIVE_TESTS=1 python -m unittest discover -s f/tests -p 'test_union*.py' -v

No HPC data is accessed. An isolated temporary directory is removed on success.
The fixture validates classifications, pre-confidence union, hit-group handling,
paired LCA, per-minimizer batch/scalar oracle, and fail-closed metadata checks.
"""
from __future__ import annotations

import os
import copy
import csv
import gzip
import io
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

F_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(F_ROOT))
from native_union.install import install


def fasta(name, sequence):
    return f">{name}\n{sequence}\n"


def fastq(name, sequence):
    return f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"


@unittest.skipUnless(os.environ.get("METAG_RUN_NATIVE_TESTS") == "1",
                     "set METAG_RUN_NATIVE_TESTS=1 for native compilation tests")
class NativeUnionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="metag-native-union-")
        cls.base = Path(cls.temporary.name)
        receipt = install({"tools_dir": str(cls.base / "tools"), "setup_cores": 4})
        cls.union = Path(receipt["executable"])
        cls.src = cls.base / "tools/metag_union/source/src"
        # Production metag_build performs the same padding initialization. This
        # local oracle builder otherwise remains the original deterministic build.
        builder = cls.src / "build_db.cc"
        text = builder.read_text()
        assert text.count("  IndexOptions index_opts;") == 1
        builder.write_text(text.replace("  IndexOptions index_opts;", "  IndexOptions index_opts = {};"))
        cls.compile_log = cls.base / "fixture-compile.log"
        with cls.compile_log.open("w") as log:
            subprocess.run(["make", "-j4", "build_db", "classify"], cwd=cls.src, check=True,
                           stdout=log, stderr=subprocess.STDOUT)
        cls.original = cls.src / "classify"
        cls.taxonomy = cls.base / "ncbi"
        cls.taxonomy.mkdir()
        nodes = [(1, 1, "no rank"), (2, 1, "superkingdom"), (3, 2, "genus"),
                 (10, 3, "species"), (20, 3, "species")]
        (cls.taxonomy / "nodes.dmp").write_text("".join(
            f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            for taxid, parent, rank in nodes))
        (cls.taxonomy / "names.dmp").write_text("".join(
            f"{taxid}\t|\tFixture {taxid}\t|\t\t|\tscientific name\t|\n"
            for taxid, _, _ in nodes))
        (cls.taxonomy / "merged.dmp").write_text("")
        cls.mapping = cls.base / "global.seqid2taxid.map"
        cls.mapping.write_text("r10 10\ns10 10\nr20 20\ns20 20\n")
        rng = random.Random(135791)
        cls.unique10 = "".join(rng.choices("ACGT", k=6000))
        cls.unique20 = "".join(rng.choices("ACGT", k=6000))
        cls.shared = "".join(rng.choices("ACGT", k=2500))
        fa = fasta("r10", cls.unique10) + fasta("s10", cls.shared)
        fb = fasta("r20", cls.unique20) + fasta("s20", cls.shared)
        cls.a, cls.b, cls.mono = (cls.base / name for name in ("a", "b", "monolithic"))
        for db, data in ((cls.a, fa), (cls.b, fb), (cls.mono, fa + fb)):
            db.mkdir()
            command = [str(cls.src / "build_db"), "-H", str(db / "hash.k2d"),
                       "-m", str(cls.mapping), "-t", str(db / "taxo.k2d"),
                       "-n", str(cls.taxonomy), "-o", str(db / "opts.k2d"),
                       "-k", "35", "-l", "31", "-c", "200000", "-p", "2"]
            subprocess.run(command, input=data, text=True, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cls.hashlist = cls.base / "hashes.txt"
        cls.write_list(cls.hashlist, (cls.a, cls.b))
        assert (cls.a / "taxo.k2d").read_bytes() == (cls.b / "taxo.k2d").read_bytes()
        assert (cls.a / "opts.k2d").read_bytes() == (cls.b / "opts.k2d").read_bytes()

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    @staticmethod
    def write_list(path, databases):
        path.write_text("".join(str(db / "hash.k2d") + "\n" for db in databases))

    def run_query(self, data, *, baseline=False, paired=False, hashes=None,
                  confidence="0.1", groups="2", extra=()):
        db = self.mono if baseline else self.a
        executable = self.original if baseline else self.union
        cmd = [str(executable), "-H", str(db / "hash.k2d"), "-t", str(db / "taxo.k2d"),
               "-o", str(db / "opts.k2d"), "-M", "-p", "2", "-T", confidence,
               "-g", groups, "-Q", "0"]
        if not baseline:
            cmd += ["--metag-hashes", str(hashes or self.hashlist)]
        if paired:
            cmd += ["-S"]
        return subprocess.run(cmd + list(extra), input=data, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def records(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [row.split("\t") for row in result.stdout.splitlines()]
        self.assertTrue(all(len(row) == 5 for row in rows), result.stdout)
        return {row[1]: row for row in rows}

    def test_single_end_matches_fixture_monolith_and_rejects_weak_hits(self):
        reads = [self.unique10[1000:1300], self.unique20[1000:1300], self.shared[500:800],
                 "N" * 150, self.unique10[:150] + "N" * 35 + self.unique20[:150],
                 self.unique10[900:940] + "N" * 1460, self.unique10[1600:1635]]
        data = "".join(fastq(str(i), sequence) for i, sequence in enumerate(reads))
        union = self.records(self.run_query(data))
        mono = self.records(self.run_query(data, baseline=True))
        self.assertEqual(union, mono)
        # The chimeric read's original minimizer counts are slightly asymmetric
        # at the ambiguous-run boundary; native scoring chooses species 20.
        # Its whole five-column output must still equal the monolithic control.
        self.assertEqual([int(union[str(i)][2]) for i in range(len(reads))],
                         [10, 20, 3, 0, 20, 0, 0])
        relaxed = self.records(self.run_query(data, confidence="0", groups="1"))
        self.assertEqual(relaxed["5"][2], "10")
        self.assertEqual(relaxed["6"][2], "10")

    def test_pairs_combine_evidence_once_and_preserve_global_hit_groups(self):
        pairs = [(self.unique10[100:250], self.unique10[800:950]),
                 (self.unique10[100:250], self.unique20[800:950]),
                 (self.shared[100:250], self.shared[800:950]),
                 ("N" * 150, "N" * 150),
                 (self.unique10[100:135], "N" * 150),
                 (self.unique10[100:135], self.unique10[800:835])]
        data = "".join(fastq(f"{i}/1", a) + fastq(f"{i}/2", b)
                       for i, (a, b) in enumerate(pairs))
        union = self.records(self.run_query(data, paired=True))
        mono = self.records(self.run_query(data, paired=True, baseline=True))
        self.assertEqual(union, mono)
        self.assertEqual([int(union[str(i)][2]) for i in range(len(pairs))], [10, 3, 3, 0, 0, 10])
        self.assertTrue(all("|" in row[3] and "|:|" in row[4] for row in union.values()))

    def test_order_of_categories_does_not_change_evidence(self):
        reversed_list = self.base / "reversed.txt"
        self.write_list(reversed_list, (self.b, self.a))
        data = fastq("1", self.shared[100:400]) + fastq("2", self.unique10[100:400])
        self.assertEqual(self.records(self.run_query(data)),
                         self.records(self.run_query(data, hashes=reversed_list)))

    def test_rejects_taxonomy_options_truncation_duplicate_and_empty_lists(self):
        data = fastq("1", self.unique10[100:400])
        for target in ("taxo.k2d", "opts.k2d", "hash.k2d"):
            bad = self.base / ("bad-" + target)
            shutil.copytree(self.b, bad)
            path = bad / target
            content = bytearray(path.read_bytes())
            if target == "hash.k2d":
                del content[-4:]
            else:
                content[-1] ^= 1
            path.write_bytes(content)
            listed = self.base / (target + ".txt")
            self.write_list(listed, (self.a, bad))
            result = self.run_query(data, hashes=listed)
            self.assertNotEqual(result.returncode, 0, target)
            self.assertEqual(result.stdout, "", target)
        for name, content in (("empty", ""), ("duplicate", str(self.a / "hash.k2d") + "\n" +
                                              str(self.a / "hash.k2d") + "\n")):
            listed = self.base / (name + ".txt")
            listed.write_text(content)
            result = self.run_query(data, hashes=listed)
            self.assertNotEqual(result.returncode, 0, name)
            self.assertEqual(result.stdout, "", name)

    def published_catalog(self, name):
        import partitioned_build
        from hpc_common import atomic_json, fingerprint

        root = self.base / name
        root.mkdir()
        cfg = {"database_layout": "partitioned", "classifier_engine": "kraken",
               "library_dir": str(root / "union_library"),
               "tools_dir": str(self.base / "tools"), "reference_raw": str(root / "references"),
               "metadata_dir": str(self.taxonomy), "reads_raw": str(root / "reads"),
               "work_dir": str(root / "work"), "output_dir": str(root / "outputs"),
               "classify_threads": 2, "minimum_hit_groups": 2, "confidence": 0.1,
               "minimum_base_quality": 0, "min_free_gib": 0, "output_format": "tar",
               "classify_cores": 8, "classify_workers": 1, "classify_memory_gib": 2,
               "system_reserve_gib": 0.1, "worker_memory_gib": 0.1, "memory_accounting": "shared"}
        for key in ("reference_raw", "reads_raw", "work_dir"):
            Path(cfg[key]).mkdir()
        index = partitioned_build.index_root(cfg)
        index.mkdir(parents=True)
        preparation = {"identity": fingerprint({"fixture": "native-union"}),
                       "groups": {"microbes": {}, "vertebrates": {}},
                       "taxonomy_dir": str(self.taxonomy)}
        for category, source in (("microbes", self.a), ("vertebrates", self.b)):
            folder = index / category
            shutil.copytree(source, folder)
            files = {name: partitioned_build.record(folder / name)
                     for name in ("hash.k2d", "taxo.k2d", "opts.k2d")}
            context = {"category": category, "db_dir": str(folder),
                       "build_identity": fingerprint({"fixture_category": category}),
                       "preparation_identity": preparation["identity"],
                       "files": files, "resource": {"capacity": 200000},
                       "parameters": dict(partitioned_build.PARAMETERS)}
            atomic_json(folder / "database.json", {"schema": partitioned_build.SCHEMA,
                        "status": "complete", "identity": fingerprint(context), "context": context})
        with patch.object(partitioned_build, "_preparation", return_value=preparation):
            snapshot = partitioned_build.publish(cfg)
            self.assertEqual(partitioned_build.publish(cfg), snapshot)
        return cfg, snapshot

    def test_published_catalog_ordinal_split_conserves_1001_pair_bytes_and_resume(self):
        import hpc_query
        from hpc_common import effective_workers

        cfg, snapshot = self.published_catalog("end-to-end")
        self.assertEqual(set(snapshot["context"]["absent_categories"]),
                         {"viruses", "plants", "invertebrates"})
        self.assertEqual(effective_workers(cfg, snapshot)["classify_workers"], 1)
        records, truth = [[], []], []
        for ordinal in range(1001):
            kind = ordinal % 5
            truth.append((10, 20, 3, 0, 3)[kind])
            sequences = ((self.unique10, self.unique10), (self.unique20, self.unique20),
                         (self.shared, self.shared), ("N" * 1000, "N" * 1000),
                         (self.unique10, self.unique20))[kind]
            for mate in range(2):
                seq = sequences[mate][100 + mate * 200:250 + mate * 200].encode()
                raw = (b"@duplicate/" + str(mate + 1).encode() + b" repeated id\r\n" + seq +
                       b"\r\n+keep this original comment\r\n" + b"I" * 150 +
                       (b"" if ordinal == 1000 else b"\r\n"))
                records[mate].append(raw)
        paths = [Path(cfg["reads_raw"]) / f"R{mate + 1}.fastq" for mate in range(2)]
        for path, mate_records in zip(paths, records):
            path.write_bytes(b"".join(mate_records))
        sample = {"dataset": "fixture", "sample": "sample", "sample_id": "fixture/sample",
                  "r1": str(paths[0]), "r2": str(paths[1])}
        result = hpc_query.run_sample(cfg, sample, snapshot)
        self.assertEqual(result["status"], "complete")
        with tarfile.open(Path(result["output"]) / "species.tar") as archive:
            with gzip.GzipFile(fileobj=archive.extractfile("read-assignments.tsv.gz")) as stream:
                rows = list(csv.DictReader(io.TextIOWrapper(stream), delimiter="\t"))
            self.assertEqual([int(row["assigned_taxid"]) for row in rows], truth)
            for tid in (0, 3, 10, 20):
                bin_path = next(row["bin"] for row in rows if int(row["assigned_taxid"]) == tid)
                for mate in range(2):
                    expected = b"".join(records[mate][i] for i, taxid in enumerate(truth) if taxid == tid)
                    with gzip.GzipFile(fileobj=archive.extractfile(f"{bin_path}/PE.R{mate + 1}.fastq.gz")) as stream:
                        self.assertEqual(stream.read(), expected)
        self.assertEqual(hpc_query.run_sample(cfg, sample, snapshot)["status"], "skipped")
        # A syntactically valid same-size mutation must not reuse prior results.
        original = paths[0].read_bytes()
        changed = original.replace(b"repeated id", b"replaced id", 1)
        self.assertEqual(len(original), len(changed))
        paths[0].write_bytes(changed)
        with self.assertRaisesRegex(RuntimeError, "identity conflicts"):
            hpc_query.run_sample(cfg, sample, snapshot)

    def test_catalog_rejects_self_consistent_incomplete_or_duplicate_category_lists(self):
        import partitioned_build
        from hpc_common import atomic_json, fingerprint

        cfg, snapshot = self.published_catalog("malformed-catalog")
        db = partitioned_build.index_root(cfg)
        for defect in ("missing_hash", "duplicate_hash", "missing_category_file", "missing_absent_record"):
            with self.subTest(defect=defect):
                bad = copy.deepcopy(snapshot)
                context = bad["context"]
                if defect == "missing_hash":
                    context["hashes"] = context["hashes"][:1]
                elif defect == "duplicate_hash":
                    context["hashes"].append(context["hashes"][0])
                elif defect == "missing_category_file":
                    del context["files"]["vertebrates/opts.k2d"]
                else:
                    context["absent_categories"].pop()
                for key in ("layout", "db_dir", "taxonomy_dir", "hashes", "files", "taxonomy_files", "parameters"):
                    bad[key] = context[key]
                bad["identity"] = fingerprint(context)
                atomic_json(db / "database.json", bad)
                # Recomputing the catalog identity must not make an incomplete
                # scientific reference set acceptable to the independent reader.
                with self.assertRaisesRegex(RuntimeError, "every recorded reference category|exact category|every absent category"):
                    partitioned_build.validate_snapshot(bad)
        atomic_json(db / "database.json", snapshot)
        partitioned_build.validate_snapshot(snapshot)
        missing = db / "vertebrates/hash.k2d"
        missing.rename(missing.with_suffix(".absent"))
        with self.assertRaises((OSError, RuntimeError)):
            partitioned_build.validate_snapshot(snapshot)

    def test_runtime_touch_requires_full_audit_then_reuses_verified_cache(self):
        import partitioned_build
        from hpc_common import sha256

        cfg, snapshot = self.published_catalog("verification-cache")
        db = partitioned_build.index_root(cfg)
        path = db / "microbes/hash.k2d"
        before = path.stat()
        os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000))
        with patch.object(partitioned_build, "sha256", side_effect=AssertionError("Quick validation rehashed an index")):
            with self.assertRaisesRegex(RuntimeError, "audit --full"):
                partitioned_build.validate_snapshot(snapshot)
        with patch.object(partitioned_build, "sha256", wraps=sha256) as digests:
            partitioned_build.validate_snapshot(snapshot, full=True)
        self.assertIn(path, [Path(call.args[0]) for call in digests.call_args_list])
        verified = json.loads((db / "verified-runtime.json").read_text())
        self.assertEqual(verified["identity"], snapshot["identity"])
        self.assertEqual(verified["files"][str(path)], {
            "state": partitioned_build.runtime_state(path),
            "sha256": snapshot["files"]["microbes/hash.k2d"]["sha256"]})
        with patch.object(partitioned_build, "sha256", side_effect=AssertionError("Verified state was rehashed")):
            self.assertEqual(partitioned_build.validate_snapshot(snapshot), snapshot)
            self.assertEqual(partitioned_build.validate_snapshot(snapshot), snapshot)
        self.assertEqual(json.loads((db / "database.json").read_text()), snapshot)
        # A cache for a prior state cannot conceal an equal-length data change.
        with path.open("r+b") as handle:
            handle.seek(32)
            original = handle.read(1)
            handle.seek(32)
            handle.write(bytes([original[0] ^ 1]))
        with self.assertRaisesRegex(RuntimeError, "audit --full"):
            partitioned_build.validate_snapshot(snapshot)
        with self.assertRaisesRegex(RuntimeError, "SHA256 mismatch"):
            partitioned_build.validate_snapshot(snapshot, full=True)

    def test_batched_lookup_matches_independent_scalar_minimizer_lca_oracle(self):
        harness = self.base / "oracle.cc"
        harness.write_text(r'''
#include <err.h>
#include "metag_union.h"
#include "mmscanner.h"
#include <iostream>
int main(int argc, char **argv) {
  using namespace kraken2;
  std::string a = std::string(argv[1]) + "/hash.k2d", b = std::string(argv[2]) + "/hash.k2d";
  Taxonomy taxonomy(std::string(argv[1]) + "/taxo.k2d", true);
  CompactHashTable<CompactHashCell> ah(a, true), bh(b, true);
  metag_union::UnionHashTable combined({a, b}, taxonomy);
  IndexOptions options{};
  std::ifstream opt(std::string(argv[1]) + "/opts.k2d", std::ios::binary);
  opt.read(reinterpret_cast<char *>(&options), sizeof(options));
  MinimizerScanner scanner(options.k, options.l, options.spaced_seed_mask,
                           true, options.toggle_mask, options.revcom_version);
  std::vector<hkey_t> keys;
  std::string sequence;
  while (std::getline(std::cin, sequence)) {
    scanner.LoadSequence(sequence);
    uint64_t *key;
    while ((key = scanner.NextMinimizer())) if (!scanner.is_ambiguous()) keys.push_back(*key);
  }
  for (uint64_t i = 0; i < 1000; ++i) keys.push_back(MurmurHash3(i) & ((uint64_t(1) << 60) - 1));
  std::vector<hvalue_t> actual(keys.size());
  combined.GetBatch(keys.data(), actual.data(), keys.size());
  size_t shared = 0, absent = 0;
  for (size_t i = 0; i < keys.size(); ++i) {
    hvalue_t x = ah.Get(keys[i]), y = bh.Get(keys[i]);
    hvalue_t expected = static_cast<hvalue_t>(taxonomy.LowestCommonAncestor(x, y));
    if (actual[i] != expected || combined.Get(keys[i]) != expected) return 1;
    if (x && y && x != y && taxonomy.nodes()[expected].external_id == 3) shared++;
    if (!expected) absent++;
  }
  if (!shared || !absent || keys.size() < 1000) return 2;
  std::cout << keys.size() << " " << shared << " " << absent << "\n";
}
''')
        executable = self.base / "oracle"
        subprocess.run(["g++", "-std=c++11", "-O2", "-fopenmp", "-DLINEAR_PROBING",
                        "-I", str(self.src), str(harness),
                        *[str(self.src / (name + ".o")) for name in
                          ("taxonomy", "mmap_file", "compact_hash", "mmscanner", "omp_hack")],
                        "-o", str(executable)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = subprocess.run([str(executable), str(self.a), str(self.b)],
                                input="\n".join((self.unique10, self.unique20, self.shared)) + "\n",
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr)
        keys, shared, absent = map(int, result.stdout.split())
        self.assertGreater(keys, 10000)
        self.assertGreater(shared, 100)
        self.assertGreater(absent, 100)


if __name__ == "__main__":
    unittest.main()
