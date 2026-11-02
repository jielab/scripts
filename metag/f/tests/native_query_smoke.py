#!/usr/bin/env python3
"""Check all three real native adapters against a prebuilt two-species fixture.

This is an explicit integration check, not part of unittest discovery.  The
fixture config must identify a complete Kraken index, its Kun-peng conversion,
and installed native tools.  Its reference_raw contains species10.fna and
species20.fna with at least 100,000 bases each, mapped to taxids 10 and 20.

Example:
  python3 native_query_smoke.py --code-root /path/to/metag \
      --fixture-config /path/to/fixture/config.json

All new reads and outputs go to a fresh temporary directory.  The fixture
references, indexes, configuration and native executables are read-only inputs.
"""
import argparse
import collections
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--fixture-config", type=Path, required=True)
    parser.add_argument("--output-parent", type=Path)
    parser.add_argument("--fragments", type=int, default=1001)
    args = parser.parse_args()
    if args.fragments < 3:
        parser.error("--fragments must be at least 3")
    sys.path.insert(0, str(args.code_root.resolve() / "f"))
    from hpc_common import load_config
    import hpc_query

    root = Path(tempfile.mkdtemp(prefix="metag_native_smoke_", dir=args.output_parent))
    cfg_base = load_config(args.fixture_config)
    cfg_base.update(reads_raw=str(root / "reads"), work_dir=str(root / "work"),
                    min_free_gib=0, classify_threads=4, output_format="tar")
    Path(cfg_base["reads_raw"]).mkdir()
    refs = {}
    for taxid in (10, 20):
        path = Path(cfg_base["reference_raw"]) / f"species{taxid}.fna"
        refs[taxid] = b"".join(path.read_bytes().splitlines()[1:])
        if len(refs[taxid]) < 100000:
            raise RuntimeError(f"Expected the two-species test fixture, short reference: {path}")

    records = [[], []]
    truth = []
    for ordinal in range(args.fragments):
        taxid = (10, 20, 0)[ordinal % 3]
        truth.append(taxid)
        for mate in range(2):
            offset = (ordinal * 61 + mate * 167) % 99000
            seq = refs[taxid][offset:offset + 150] if taxid else b"N" * 150
            record = (b"@duplicate/" + str(mate + 1).encode() + b" same name repeats\r\n" +
                      seq + b"\r\n+original comment\r\n" + b"I" * 150 +
                      (b"" if ordinal == args.fragments - 1 else b"\r\n"))
            records[mate].append(record)
    for mate in range(2):
        path = Path(cfg_base["reads_raw"]) / f"sample_R{mate + 1}.fastq.gz"
        with gzip.open(path, "wb", compresslevel=1) as handle:
            handle.write(b"".join(records[mate]))
    sample = {"sample_id": "fixture/sample", "dataset": "fixture", "sample": "sample",
              "r1": str(Path(cfg_base["reads_raw"]) / "sample_R1.fastq.gz"),
              "r2": str(Path(cfg_base["reads_raw"]) / "sample_R2.fastq.gz")}

    results = []
    for engine in ("kraken", "kunpeng_full", "kunpeng_memory"):
        cfg = dict(cfg_base, classifier_engine=engine, output_dir=str(root / "outputs" / engine))
        started = time.monotonic()
        result = hpc_query.run_sample(cfg, sample)
        elapsed = time.monotonic() - started
        archive = Path(result["output"]) / "species.tar"
        with tarfile.open(archive, "r") as tar:
            with gzip.GzipFile(fileobj=tar.extractfile("read-assignments.tsv.gz")) as handle:
                rows = list(csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8"), delimiter="\t"))
            assignments = [int(row["species_taxid"]) for row in rows]
            if assignments != truth:
                raise AssertionError(f"{engine}: native assignments differ from fixture truth")
            for taxid in (0, 10, 20):
                bin_path = next(row["bin"] for row in rows if int(row["species_taxid"]) == taxid)
                for mate in range(2):
                    with gzip.GzipFile(fileobj=tar.extractfile(f"{bin_path}/PE.R{mate + 1}.fastq.gz")) as handle:
                        expected = b"".join(records[mate][i] for i, tid in enumerate(truth) if tid == taxid)
                        if handle.read() != expected:
                            raise AssertionError(f"{engine}: original FASTQ bytes changed in {bin_path}")
        repeated = hpc_query.run_sample(cfg, sample)
        if repeated["status"] != "skipped":
            raise AssertionError(f"{engine}: completed result was not reused")
        # A valid, equal-length replacement must never reuse the first species
        # result. Exercise this with real native engines, not mocked stdout.
        resume_path = Path(cfg_base["reads_raw"]) / f"{engine}.resume.fastq"
        def one_read(taxid):
            return b"@same-id\n" + refs[taxid][:150] + b"\n+\n" + b"I" * 150 + b"\n"
        original_bytes, replacement_bytes = one_read(10), one_read(20)
        if len(original_bytes) != len(replacement_bytes):
            raise AssertionError("Invalid same-size test fixture")
        resume_path.write_bytes(original_bytes)
        resume_sample = dict(sample_id="fixture/resume", dataset="fixture", sample="resume",
                             r1=str(resume_path), r2=None)
        first = hpc_query.run_sample(cfg, resume_sample)
        def assigned_species(published):
            with (Path(published["output"]) / "counts.tsv").open() as handle:
                return [int(row["species_taxid"]) for row in csv.DictReader(handle, delimiter="\t")]
        if assigned_species(first) != [10]:
            raise AssertionError(f"{engine}: first single-read truth failed")
        resume_path.write_bytes(replacement_bytes)
        try:
            hpc_query.run_sample(cfg, resume_sample)
        except RuntimeError as error:
            if "identity conflicts" not in str(error):
                raise
        else:
            raise AssertionError(f"{engine}: silently reused a same-size changed FASTQ")
        fresh_cfg = dict(cfg, output_dir=str(root / "fresh-outputs" / engine))
        fresh = hpc_query.run_sample(fresh_cfg, resume_sample)
        if assigned_species(fresh) != [20]:
            raise AssertionError(f"{engine}: replacement read truth failed in fresh output")
        verification = {"engine": engine, "fragments": len(rows),
                        "truth_counts": dict(collections.Counter(truth)),
                        "exact_bin_fastq_bytes": True, "duplicate_ids_preserved": True,
                        "CRLF_plus_comments_final_newline_preserved": True,
                        "same_size_valid_input_replacement_rejected": True,
                        "fresh_replacement_species": 20,
                        "resume": repeated["status"], "seconds": elapsed}
        print(json.dumps(verification, sort_keys=True), flush=True)
        results.append(dict(verification, summary=result["summary"]))
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(f"RESULT_ROOT {root}")


if __name__ == "__main__":
    main()
