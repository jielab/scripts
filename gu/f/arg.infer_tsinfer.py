#!/usr/bin/env python3
"""Build a dated tsinfer/tsdate tree sequence for TRACE from a filtered phased VCF.

The workflow uses explicit generate_ancestors -> match_ancestors -> match_samples
stages rather than tsinfer.infer(). This makes memory use and failure points easier to
inspect, and it allows mmap-backed ancestor generation on local Linux storage.

The input should be produced by ``f/arg.prep_vcf.sh``. It must contain only
biallelic variants, phased GTs, and a normalized INFO/AA field. INFO/AA, rather
than VCF REF, is used as the ancestral state.

Large VCF-Zarr intermediates are created below a Linux temporary directory
(default: /tmp), not beside the final output on /mnt/d. This is important in
WSL because Zarr contains many small files and is extremely slow to create or
remove on a Windows-mounted filesystem.
"""

from __future__ import annotations

import argparse
from collections import Counter
import inspect
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import time
import warnings
from typing import Any


def die(message: str) -> "NoReturn":
    raise SystemExit(message)


def command_path(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        die(f"Required command not found: {name}")
    return path


def version_of(module: Any) -> str:
    return str(getattr(module, "__version__", "unknown"))


def extract_sample_ids(vcf_path: Path) -> list[str]:
    """Read sample IDs with bcftools, avoiding Python gzip/header edge cases."""
    result = subprocess.run(
        [command_path("bcftools"), "query", "-l", str(vcf_path)],
        check=True,
        text=True,
        capture_output=True,
    )
    samples = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not samples:
        die(f"No samples found in VCF: {vcf_path}")
    if len(samples) != len(set(samples)):
        die(f"Duplicate sample IDs found in VCF: {vcf_path}")
    return samples


def check_phasing(vcf_path: Path, max_records: int = 2000) -> dict[str, int]:
    """Check diploid heterozygotes; haploid calls are phase-complete by definition."""
    import pysam

    records = 0
    called = 0
    heterozygous = 0
    unphased_heterozygous = 0
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            records += 1
            for call in record.samples.values():
                gt = call.get("GT")
                if gt is None or any(a is None for a in gt):
                    continue
                called += 1
                if len(gt) < 2:
                    continue
                if len(set(gt)) > 1:
                    heterozygous += 1
                    if not call.phased:
                        unphased_heterozygous += 1
                        if unphased_heterozygous <= 5:
                            print(
                                f"UNPHASED: {record.contig}:{record.pos} GT={gt}",
                                flush=True,
                            )
            if records >= max_records:
                break
    if records == 0:
        die(f"VCF has no records: {vcf_path}")
    if unphased_heterozygous:
        die(
            f"VCF is not fully phased: found {unphased_heterozygous} unphased "
            f"heterozygous calls in the first {records} records"
        )
    return {
        "records_checked": records,
        "called_genotypes_checked": called,
        "heterozygous_genotypes_checked": heterozygous,
        "unphased_heterozygous_checked": unphased_heterozygous,
    }


def scalar_string(value: Any) -> str:
    """Convert Zarr string/byte/scalar/length-one arrays to a plain string."""
    if hasattr(value, "tolist"):
        value = value.tolist()
    while isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value)


def normalize_aa(value: Any) -> str:
    aa = scalar_string(value).strip().upper()
    aa = aa.split("|", 1)[0].split(",", 1)[0]
    return aa


def build_ancestral_array(zarr_path: Path, chunk_size: int = 250_000) -> tuple[str, int]:
    """Create a validated one-base ancestral-state array from variant_AA."""
    import numpy as np
    import zarr

    root = zarr.open(str(zarr_path), mode="r+")
    required = ("variant_AA", "variant_allele")
    missing = [name for name in required if name not in root]
    if missing:
        die(
            "VCF-Zarr is missing "
            + ", ".join(missing)
            + ". The VCF must retain INFO/AA; REF is not an acceptable fallback."
        )
    aa_source = root["variant_AA"]
    alleles_source = root["variant_allele"]
    n_variants = int(aa_source.shape[0])
    if n_variants == 0:
        die("VCF-Zarr contains no variants")
    out_name = "variant_ancestral_state"
    if out_name in root:
        del root[out_name]
    chunks = getattr(aa_source, "chunks", None)
    out_chunk = int(chunks[0]) if chunks and chunks[0] else min(chunk_size, n_variants)
    out = root.create_dataset(
        out_name,
        shape=(n_variants,),
        chunks=(max(1, out_chunk),),
        dtype="<U1",
        overwrite=True,
    )

    invalid: list[str] = []
    for start in range(0, n_variants, chunk_size):
        end = min(start + chunk_size, n_variants)
        aa_block = aa_source[start:end]
        allele_block = alleles_source[start:end]
        normalized = np.empty(end - start, dtype="<U1")
        for offset, raw_aa in enumerate(aa_block):
            aa = normalize_aa(raw_aa)
            alleles = {
                scalar_string(x).strip().upper()
                for x in allele_block[offset]
                if scalar_string(x).strip() not in ("", ".")
            }
            if aa not in {"A", "C", "G", "T"} or aa not in alleles:
                if len(invalid) < 10:
                    invalid.append(
                        f"variant_index={start + offset} AA={aa!r} alleles={sorted(alleles)}"
                    )
                normalized[offset] = "N"
            else:
                normalized[offset] = aa
        out[start:end] = normalized
    if invalid:
        die(
            "Invalid ancestral states remain after VCF preparation. Examples: "
            + "; ".join(invalid)
        )
    return out_name, n_variants


def vcf_zarr_ploidy(zarr_path: Path) -> int:
    """Read the fixed GT ploidy dimension written by vcf2zarr."""
    import zarr

    root = zarr.open(str(zarr_path), mode="r")
    if "call_genotype" not in root:
        die("VCF-Zarr is missing call_genotype")
    shape = root["call_genotype"].shape
    if len(shape) != 3 or int(shape[2]) < 1:
        die(f"Unexpected call_genotype shape in VCF-Zarr: {shape}")
    return int(shape[2])


def metadata_dict(obj: Any) -> dict[str, Any]:
    metadata = getattr(obj, "metadata", None)
    if isinstance(metadata, dict):
        return metadata
    if isinstance(metadata, bytes):
        try:
            return json.loads(metadata.decode())
        except Exception:
            return {}
    return {}


def write_sample_map(tree_sequence: Any, sample_ids: list[str], output_path: Path) -> None:
    """Write actual tree node IDs and map haplotypes to VCF samples."""
    seen: dict[str, int] = {}
    sample_id_set = set(sample_ids)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("tree_node_id\tsample\thaplotype\n")
        for node_id in tree_sequence.samples():
            node = tree_sequence.node(int(node_id))
            candidates: list[str] = []
            objects = [node]
            if node.individual >= 0:
                objects.append(tree_sequence.individual(node.individual))
            for obj in objects:
                metadata = metadata_dict(obj)
                for key in (
                    "variant_data_sample_id",
                    "sample",
                    "sample_id",
                    "individual",
                    "individual_id",
                    "name",
                    "id",
                ):
                    if key in metadata:
                        candidates.append(str(metadata[key]))
            sample = next((x for x in candidates if x in sample_id_set), None)
            if sample is None and 0 <= node.individual < len(sample_ids):
                sample = sample_ids[node.individual]
            if sample is None:
                sample = candidates[0] if candidates else f"node_{node_id}"
            seen[sample] = seen.get(sample, 0) + 1
            handle.write(f"{int(node_id)}\t{sample}\t{seen[sample]}\n")


def dump_atomic(tree_sequence: Any, output: Path) -> None:
    temp_output = output.with_name(f"{output.name}.part.{os.getpid()}")
    try:
        tree_sequence.dump(str(temp_output))
        os.replace(temp_output, output)
    finally:
        temp_output.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Infer and date an ARG for TRACE from a phased INFO/AA VCF"
    )
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument("--qc-json", type=Path)
    parser.add_argument("--chr", default="")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--tmp-root", type=Path, default=Path("/tmp"))
    parser.add_argument("--recombination-rate", type=float, default=1.0e-8)
    parser.add_argument("--mutation-rate", type=float, default=1.25e-8)
    parser.add_argument("--population-size", type=float, default=30_000)
    parser.add_argument(
        "--tsdate-method",
        default="variational_gamma",
        choices=("variational_gamma", "inside_outside", "maximization"),
    )
    parser.add_argument("--mismatch-ratio", type=float)
    parser.add_argument("--phase-check-records", type=int, default=2000)
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument("--one-bit", action="store_true",
                        help="use 1-bit biallelic genotype encoding; requires no missing genotypes")
    parser.add_argument(
        "--sample-match-min-work",
        type=int,
        default=50_000_000,
        help=(
            "genotypes per resumable match-samples partition; 0 disables batching "
            "(default: 50000000)"
        ),
    )
    parser.add_argument("--keep-work", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        die("--threads must be positive")
    if args.sample_match_min_work < 0:
        die("--sample-match-min-work must be nonnegative")
    for name, value in (
        ("recombination rate", args.recombination_rate),
        ("mutation rate", args.mutation_rate),
        ("population size", args.population_size),
    ):
        if value <= 0:
            die(f"{name} must be positive")
    args.vcf = args.vcf.resolve()
    if not args.vcf.is_file() or args.vcf.stat().st_size == 0:
        die(f"Input VCF is missing or empty: {args.vcf}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    tmp_filesystem = subprocess.check_output(
        [command_path("stat"), "-f", "-c", "%T", str(args.tmp_root)], text=True
    ).strip()
    if tmp_filesystem in {"v9fs", "9p", "drvfs"}:
        die(
            "ARG --tmp-root must be on the local WSL Linux filesystem, not "
            f"{tmp_filesystem}: {args.tmp_root}. LMDB requires a 1 TiB sparse map; "
            "use /tmp or set GU_ARG_TMP_ROOT to an ext4 path."
        )
    sample_map = args.sample_map or args.out.with_suffix(".sample_map.tsv")
    qc_json = args.qc_json or args.out.with_suffix(".arg_qc.json")
    sample_map.parent.mkdir(parents=True, exist_ok=True)
    qc_json.parent.mkdir(parents=True, exist_ok=True)

    warnings.filterwarnings(
        "ignore",
        message=r"The LMDBStore is deprecated and will be removed in a Zarr-Python version 3.*",
        category=FutureWarning,
    )
    import tsdate
    import tsinfer
    import tskit
    import zarr

    sample_ids = extract_sample_ids(args.vcf)
    phase_qc = check_phasing(args.vcf, args.phase_check_records)
    # Checkpoint identity follows the requested output, not input metadata.
    work = args.tmp_root / f"gu_arg_{args.out.stem}"
    work.mkdir(parents=True, exist_ok=True)
    started = time.time()
    completed = False
    try:
        vcz_path = work / f"{args.out.stem}.vcz"
        inferred_path = work / f"{args.out.stem}.inferred.trees"
        processed_path = work / f"{args.out.stem}.preprocessed.trees"

        print(f"[ARG] Work directory: {work}", flush=True)
        if (vcz_path / ".zmetadata").is_file():
            print(f"[ARG] Reuse VCF-Zarr checkpoint: {vcz_path}", flush=True)
        else:
            if vcz_path.exists():
                shutil.rmtree(vcz_path)
            print(f"[ARG] Convert VCF to VCF-Zarr: {vcz_path}", flush=True)
            subprocess.run(
                [command_path("vcf2zarr"), "convert", str(args.vcf), str(vcz_path)],
                check=True,
            )
        input_ploidy = vcf_zarr_ploidy(vcz_path)
        x_identities = None
        if args.chr.upper().removeprefix('CHR') in ('X', '23'):
            from arg_x import haploid_zarr
            haploid_path = work / f"{args.out.stem}.haploid.vcz"
            x_identities = haploid_zarr(vcz_path, haploid_path)
            vcz_path = haploid_path
            print(f"[ARG] Non-PAR X: {len(sample_ids)} individuals, {len(x_identities)} real haplotypes", flush=True)
        ancestral_array, n_input_variants = build_ancestral_array(vcz_path)
        print(f"[ARG] VCF ploidy: {input_ploidy}", flush=True)

        data = tsinfer.VariantData(
            str(vcz_path), ancestral_state=ancestral_array
        )
        # Explicit staged inference is easier to diagnose than the one-call
        # tsinfer.infer() path and permits mmap-backed ancestor generation.
        ancestor_path = work / f"{args.out.stem}.ancestors"
        ancestor_done = work / f"{args.out.stem}.ancestors.done"
        if ancestor_path.is_file() and ancestor_done.is_file():
            print(f"[ARG] Reuse ancestor-data checkpoint: {ancestor_path}", flush=True)
            ancestor_data = tsinfer.load(str(ancestor_path))
        else:
            ancestor_path.unlink(missing_ok=True)
            ancestor_done.unlink(missing_ok=True)
            print("[ARG] tsinfer.generate_ancestors", flush=True)
            gen_kwargs: dict[str, Any] = {
                "path": str(ancestor_path),
                "num_threads": args.threads,
                "mmap_temp_dir": str(work),
            }
            if args.one_bit:
                gen_kwargs["genotype_encoding"] = tsinfer.GenotypeEncoding.ONE_BIT
            ancestor_data = tsinfer.generate_ancestors(data, **gen_kwargs)
            ancestor_done.write_text("complete\n")

        match_kwargs: dict[str, Any] = {
            "recombination_rate": args.recombination_rate,
            "num_threads": args.threads,
        }
        if args.mismatch_ratio is not None:
            match_kwargs["mismatch_ratio"] = args.mismatch_ratio
        ancestors_ts_path = work / f"{args.out.stem}.ancestors.trees"
        if ancestors_ts_path.is_file():
            print(f"[ARG] Reuse ancestor-tree checkpoint: {ancestors_ts_path}", flush=True)
            ancestors_ts = tskit.load(str(ancestors_ts_path))
        else:
            print("[ARG] tsinfer.match_ancestors", flush=True)
            ancestors_ts = tsinfer.match_ancestors(data, ancestor_data, **match_kwargs)
            dump_atomic(ancestors_ts, ancestors_ts_path)

        if inferred_path.is_file():
            print(f"[ARG] Reuse matched-samples checkpoint: {inferred_path}", flush=True)
            inferred_ts = tskit.load(str(inferred_path))
        elif args.sample_match_min_work == 0:
            print("[ARG] tsinfer.match_samples", flush=True)
            inferred_ts = tsinfer.match_samples(data, ancestors_ts, **match_kwargs)
            dump_atomic(inferred_ts, inferred_path)
        else:
            batch_dir = work / "match_samples_batch"
            metadata_path = batch_dir / "metadata.json"
            batch_kwargs: dict[str, Any] = {
                "recombination_rate": args.recombination_rate,
            }
            if args.mismatch_ratio is not None:
                batch_kwargs["mismatch_ratio"] = args.mismatch_ratio
            if not metadata_path.is_file():
                print("[ARG] tsinfer.match_samples batch init", flush=True)
                descriptor = tsinfer.match_samples_batch_init(
                    str(batch_dir),
                    str(vcz_path),
                    ancestral_array,
                    str(ancestors_ts_path),
                    args.sample_match_min_work,
                    **batch_kwargs,
                )
                num_partitions = int(descriptor.num_partitions)
            else:
                metadata = json.loads(metadata_path.read_text())
                num_partitions = int(metadata["num_partitions"])
                print(
                    f"[ARG] Reuse match-samples batch checkpoint: {batch_dir}",
                    flush=True,
                )
            for partition_index in range(num_partitions):
                partition_path = batch_dir / f"partition_{partition_index}.pkl"
                if partition_path.is_file() and partition_path.stat().st_size > 0:
                    try:
                        with partition_path.open("rb") as partition_file:
                            pickle.load(partition_file)
                        continue
                    except (EOFError, OSError, pickle.UnpicklingError):
                        partition_path.unlink(missing_ok=True)
                print(
                    "[ARG] tsinfer.match_samples partition "
                    f"{partition_index + 1}/{num_partitions}",
                    flush=True,
                )
                tsinfer.match_samples_batch_partition(
                    str(batch_dir), partition_index
                )
            print("[ARG] tsinfer.match_samples batch finalise", flush=True)
            inferred_ts = tsinfer.match_samples_batch_finalise(str(batch_dir))
            dump_atomic(inferred_ts, inferred_path)

        if x_identities is not None:
            from arg_x import restore_individuals
            inferred_ts = restore_individuals(inferred_ts, x_identities)
        if hasattr(tsdate, "preprocess_ts"):
            print("[ARG] tsdate.preprocess_ts", flush=True)
            processed_ts = tsdate.preprocess_ts(inferred_ts)
        else:
            processed_ts = inferred_ts
        dump_atomic(processed_ts, processed_path)

        date_kwargs: dict[str, Any] = {
            "mutation_rate": args.mutation_rate,
            "method": args.tsdate_method,
        }
        # Variational gamma does not accept population_size.
        if args.tsdate_method != "variational_gamma":
            date_kwargs["population_size"] = args.population_size
        date_signature = inspect.signature(tsdate.date)
        if "num_threads" in date_signature.parameters:
            date_kwargs["num_threads"] = args.threads
        if "progress" in date_signature.parameters:
            date_kwargs["progress"] = not args.no_progress
        date_parameters = (
            f"method={args.tsdate_method}, mutation_rate={args.mutation_rate:g}"
        )
        if "population_size" in date_kwargs:
            date_parameters += f", Ne={args.population_size:g}"
        print(f"[ARG] tsdate: {date_parameters}", flush=True)
        dated_ts = tsdate.date(processed_ts, **date_kwargs)
        expected_nodes = len(x_identities) if x_identities is not None else input_ploidy * len(sample_ids)
        if dated_ts.num_samples != expected_nodes:
            die(
                "ARG sample-node count does not match VCF samples × ploidy "
                f"({dated_ts.num_samples} vs {len(sample_ids)} × {input_ploidy} = {expected_nodes})"
            )
        dump_atomic(dated_ts, args.out)

        map_part = sample_map.with_name(f"{sample_map.name}.part.{os.getpid()}")
        write_sample_map(dated_ts, sample_ids, map_part)
        os.replace(map_part, sample_map)

        qc = {
            "input_vcf": str(args.vcf),
            "output_trees": str(args.out.resolve()),
            "sample_map": str(sample_map.resolve()),
            "ancestral_state_source": "INFO/AA (VCF-Zarr variant_AA)",
            "input_variants": n_input_variants,
            "vcf_samples": len(sample_ids),
            "vcf_ploidy": input_ploidy,
            "sample_ploidies": (dict(Counter(name for name, _ in x_identities))
                                if x_identities is not None else dict.fromkeys(sample_ids, input_ploidy)),
            "sample_nodes": int(dated_ts.num_samples),
            "sites": int(dated_ts.num_sites),
            "mutations": int(dated_ts.num_mutations),
            "nodes": int(dated_ts.num_nodes),
            "edges": int(dated_ts.num_edges),
            "trees": int(dated_ts.num_trees),
            "sequence_length": float(dated_ts.sequence_length),
            "parameters": {
                "inference_pipeline": "generate_ancestors/match_ancestors/match_samples",
                "one_bit": bool(args.one_bit),
                "threads": args.threads,
                "sample_match_min_work": args.sample_match_min_work,
                "recombination_rate": args.recombination_rate,
                "mutation_rate": args.mutation_rate,
                "population_size": args.population_size,
                "tsdate_method": args.tsdate_method,
                "mismatch_ratio": args.mismatch_ratio,
            },
            "phase_qc": phase_qc,
            "versions": {
                "tsinfer": version_of(tsinfer),
                "tsdate": version_of(tsdate),
                "tskit": version_of(tskit),
                "zarr": version_of(zarr),
            },
            "elapsed_seconds": round(time.time() - started, 3),
        }
        qc_part = qc_json.with_name(f"{qc_json.name}.part.{os.getpid()}")
        qc_part.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
        os.replace(qc_part, qc_json)
        print(f"[ARG] Dated ARG: {args.out}", flush=True)
        print(f"[ARG] TRACE sample map: {sample_map}", flush=True)
        print(f"[ARG] QC: {qc_json}", flush=True)
        completed = True
    finally:
        if args.keep_work:
            print(
                f"[ARG] Keeping work directory by request: {work}",
                flush=True,
            )
        else:
            shutil.rmtree(work, ignore_errors=True)
            if not completed:
                print(f"[ARG] Removed incomplete work directory: {work}", flush=True)


if __name__ == "__main__":
    main()
