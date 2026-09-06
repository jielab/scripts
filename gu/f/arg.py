#!/usr/bin/env python3
"""ARG preparation, inference, validation and TRACE export utilities."""
from __future__ import annotations
import sys

# ---- x ----
"""Non-PAR X input handling and explicit biological haplotype identities."""




def x_nonpar(position, build):
    # One-based inclusive non-PAR interval between PAR1 and PAR2.
    start, end = {'37': (2699521, 154931043), '38': (2781480, 155701382)}[str(build)]
    return start <= position <= end


def x_read_sexes(path):
    with open(path) as handle:
        header = handle.readline().lstrip('#').split()
        if 'IID' not in header or 'SEX' not in header:
            raise ValueError(f'X requires IID and SEX in PSAM: {path}')
        result = {}
        for line in handle:
            if not line.strip():
                continue
            row = dict(zip(header, line.split()))
            if row['IID'] in result or row['SEX'] not in ('1', '2'):
                raise ValueError(f'Duplicate IID or unknown X sex: {row["IID"]}')
            result[row['IID']] = int(row['SEX'])
    return result


def x_haploid_zarr(source, destination):
    """Flatten only real GT copies, never infer the -2 padding as a sample.

    Uniform haploid storage also works with tsinfer's partitioned matcher.
    Scan every GT chunk to enforce constant per-individual ploidy.
    """
    import numpy as np
    import numcodecs
    import zarr

    src = zarr.open(str(source), mode='r')
    gt = src['call_genotype']
    names = [str(x) for x in src['sample_id'][:]]
    present = np.asarray(gt[0]) != -2
    ploidies = present.sum(axis=1)
    if np.any((ploidies < 1) | (ploidies > 2)) or np.any(~present[:, 0]):
        raise ValueError('Invalid per-sample X ploidy in VCF-Zarr')
    identities = [(name, h + 1) for name, p in zip(names, ploidies) for h in range(int(p))]
    dest = zarr.open_group(str(destination), mode='w')
    for key in ('variant_position', 'variant_allele', 'variant_AA', 'variant_contig', 'contig_id', 'contig_length'):
        if key in src:
            zarr.copy(src[key], dest, name=key)
    ids = np.asarray([f'arg_hap_{i}' for i in range(len(identities))], dtype=object)
    dest.create_dataset('sample_id', data=ids, object_codec=numcodecs.VLenUTF8())
    width = len(identities)
    chunk = min(gt.chunks[0], max(1, 32_000_000 // (gt.shape[1] * gt.shape[2])))
    out = dest.create_dataset('call_genotype', shape=(gt.shape[0], width, 1),
                              chunks=(chunk, min(width, 256), 1), dtype=gt.dtype)
    for start in range(0, gt.shape[0], chunk):
        values = np.asarray(gt[start:start + chunk])
        if np.any((values != -2) != present[None, :, :]):
            raise ValueError(f'X ploidy changes within chromosome at variant chunk {start}')
        if np.any(values < -2) or np.any(values > 1):
            raise ValueError('Non-biallelic GT in prepared X VCF')
        out[start:start + len(values), :, 0] = values[:, present]
    zarr.consolidate_metadata(dest.store)
    return identities


def x_restore_individuals(ts, identities):
    """Replace storage-only haploid identities with original biological samples."""
    import numpy as np
    import tskit

    if ts.num_samples != len(identities):
        raise ValueError('X ARG sample count differs from real input haplotypes')
    tables = ts.dump_tables()
    individuals = {}
    tables.individuals.clear()
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    node_individual = np.full(ts.num_nodes, tskit.NULL, dtype=np.int32)
    for node, (name, hap) in zip(ts.samples(), identities):
        if name not in individuals:
            individuals[name] = tables.individuals.add_row(metadata={'sample': name})
        node_individual[node] = individuals[name]
    tables.nodes.individual = node_individual
    return tables.tree_sequence()

# ---- needle_x ----
"""Pack real non-PAR X copies into ASMC's two-column storage containers."""
import argparse
import csv
import json
from pathlib import Path




def needle_x_panel(keep, psam, policy, qc):
    sexes = x_read_sexes(psam)
    path = Path(keep)
    rows = path.read_text().splitlines()
    names = [r.split()[-1] for r in rows[1:] if r.strip()]
    if not set(names) <= sexes.keys():
        raise ValueError('Selected X panel has samples absent from PSAM')
    males = [name for name in names if sexes[name] == 1]
    dropped = []
    if len(males) % 2:
        if policy != 'drop-last':
            raise ValueError('ASMC requires an even number of real haplotypes. '
                             'Use --x-odd-male drop-last to exclude the last selected male, '
                             'or provide an even-male sample panel with --keep.')
        dropped = males[-1:]
        rows = [rows[0]] + [r for r in rows[1:] if r.split()[-1] not in dropped]
        path.write_text('\n'.join(rows) + '\n')
    result = dict(policy=policy, excluded_males=dropped, individuals=len(names)-len(dropped),
                  haplotypes=sum(1 if sexes[n] == 1 else 2 for n in names)-len(dropped))
    qc = Path(qc)
    text = json.dumps(result, indent=2) + '\n'
    if not qc.is_file() or qc.read_text() != text:
        temporary = qc.with_name(qc.name + '.next')
        temporary.write_text(text)
        temporary.replace(qc)
    print('Needle X panel: ' + json.dumps(result))


def needle_x_pack(prefix, psam, build, keep):
    import numpy as np
    prefix = str(prefix)
    sexes = x_read_sexes(psam)
    sample = Path(prefix + '.sample')
    names = [line.split()[1] for line in sample.read_text().splitlines()[2:] if line.strip()]
    wanted = [line.split()[-1] for line in Path(keep).read_text().splitlines() if line.strip() and not line.startswith('#')]
    if names != wanted:
        raise ValueError('X export sample order differs from selected panel')
    identities = [(name, h) for name in names for h in range(1 if sexes[name] == 1 else 2)]
    if len(identities) % 2:
        raise ValueError('Odd real X haplotype count: select an even-male panel before export')
    male = np.array([sexes[name] == 1 for name in names])
    columns = np.array([2*i+h for i, name in enumerate(names) for h in range(1 if sexes[name] == 1 else 2)])
    haps = Path(prefix + '.haps')
    tmp = Path(prefix + '.haps.x.next')
    retained = removed = 0
    try:
        with haps.open() as src, tmp.open('w') as dest:
            for line in src:
                z = line.split()
                if len(z) != 5 + 2 * len(names) or z[0].removeprefix('chr') not in ('X', '23'):
                    raise ValueError('Malformed X HAPS row')
                if not x_nonpar(int(z[2]), build):
                    removed += 1
                    continue
                # Vectorize across individuals: full-density X contains billions
                # of calls. ASCII 0/1 are 48/49; PLINK's absent male copy is '-'.
                values = np.frombuffer(''.join(z[5:]).encode('ascii'), dtype=np.uint8)
                if len(values) != 2 * len(names):
                    raise ValueError(f'Missing/unphased X HAPS at {z[2]}')
                pairs = values.reshape(-1, 2)
                a, b = pairs[:, 0], pairs[:, 1]
                valid_a = (a == 48) | (a == 49)
                valid_b = np.where(male, (b == 45) | (b == a), (b == 48) | (b == 49))
                if not np.all(valid_a & valid_b):
                    bad = int(np.flatnonzero(~(valid_a & valid_b))[0])
                    raise ValueError(f'Invalid/missing X HAPS: {names[bad]} at {z[2]}')
                gt = values[columns].tobytes().decode('ascii')
                dest.write(' '.join(['X'] + z[1:5]) + ' ' + ' '.join(gt) + '\n')
                retained += 1
        if retained < 2:
            raise ValueError('Too few non-PAR X HAPS variants')
        tmp.replace(haps)
    finally:
        tmp.unlink(missing_ok=True)
    # Pairs are file-format containers only; biological identities live in the sidecar.
    storage_ids = [f'X_pair_{i}' for i in range(len(identities)//2)]
    sample.write_text('ID_1 ID_2 missing\n0 0 0\n' + ''.join(f'{n} {n} 0\n' for n in storage_ids))
    Path(prefix + '.packed.keep').write_text('#IID\n' + '\n'.join(storage_ids) + '\n')
    with open(prefix + '.haplotypes.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('eid', 'hap'))
        writer.writerows(identities)
    print(f'Needle X packed: individuals={len(names)} real_haplotypes={len(identities)} '
          f'nonpar_variants={retained} excluded_outside_nonpar={removed}')


def needle_x_sample_map(tree, identities, output):
    import tskit


    with open(identities) as handle:
        rows = [(r['eid'], int(r['hap'])) for r in csv.DictReader(handle, delimiter='\t')]
    ts = tskit.load(tree)
    if len(rows) != ts.num_samples:
        raise ValueError('X tree sample count differs from biological haplotype map')
    restored = x_restore_individuals(ts, rows)
    temp = Path(str(tree) + '.x.next')
    restored.dump(str(temp))
    temp.replace(tree)
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(('eid', 'hap', 'sample_node'))
        writer.writerows((name, hap, int(node)) for (name, hap), node in zip(rows, restored.samples()))


def needle_x_annotate_qc(qc, identities):
    from collections import Counter
    with open(identities) as handle:
        counts = Counter(r['eid'] for r in csv.DictReader(handle, delimiter='\t'))
    path = Path(qc)
    result = json.loads(path.read_text())
    if result['haplotypes'] != sum(counts.values()):
        raise ValueError('X HAPS QC count differs from real haplotypes')
    result['storage_pairs'] = result['individuals']
    result['individuals'] = len(counts)
    result['haploid_individuals'] = sum(n == 1 for n in counts.values())
    result['diploid_individuals'] = sum(n == 2 for n in counts.values())
    result['storage'] = 'ASMC-pairs; biological identities in chrX.haplotypes.tsv'
    path.write_text(json.dumps(result, indent=2) + '\n')


def needle_x_check(tree, identities, output, psam, keep):
    import tskit
    sexes = x_read_sexes(psam)
    names = [r.split()[-1] for r in Path(keep).read_text().splitlines()[1:] if r.strip()]
    expected = [(name, h) for name in names for h in range(1 if sexes[name] == 1 else 2)]
    with open(identities) as handle:
        actual = [(r['eid'], int(r['hap'])) for r in csv.DictReader(handle, delimiter='\t')]
    if actual != expected or len(set(actual)) != len(actual):
        raise ValueError('X biological haplotype identities disagree with panel/sex')
    ts = tskit.load(tree)
    with open(output) as handle:
        mapped = [(r['eid'], int(r['hap']), int(r['sample_node'])) for r in csv.DictReader(handle, delimiter='\t')]
    if len(expected) != ts.num_samples or mapped != [(n, h, int(node)) for (n, h), node in zip(expected, ts.samples())]:
        raise ValueError('X sample map does not match tree samples and real ploidy')
    for name, hap, node in mapped:
        individual = ts.node(node).individual
        if individual < 0 or ts.individual(individual).metadata.get('sample') != name:
            raise ValueError('X tree individual identity disagrees with sample map')
    print(f'Needle X identity CHECK PASS: individuals={len(names)} real_haplotypes={len(expected)}')


def needle_x_main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('panel')
    for key in ('keep', 'psam', 'qc'):
        p.add_argument('--' + key, required=True)
    p.add_argument('--policy', choices=('error', 'drop-last'), default='drop-last')
    p = sub.add_parser('pack')
    for key in ('prefix', 'psam', 'keep', 'build'):
        p.add_argument('--' + key, required=True)
    p = sub.add_parser('map')
    for key in ('tree', 'identities', 'output'):
        p.add_argument('--' + key, required=True)
    p = sub.add_parser('check')
    for key in ('tree', 'identities', 'output', 'psam', 'keep'):
        p.add_argument('--' + key, required=True)
    p = sub.add_parser('qc')
    for key in ('qc', 'identities'):
        p.add_argument('--' + key, required=True)
    args = vars(parser.parse_args())
    command = args.pop('command')
    {'panel': needle_x_panel, 'pack': needle_x_pack, 'map': needle_x_sample_map, 'check': needle_x_check, 'qc': needle_x_annotate_qc}[command](**args)





# ---- infer ----
"""Build a dated tsinfer/tsdate tree sequence for TRACE from a filtered phased VCF.

The workflow uses explicit generate_ancestors -> match_ancestors -> match_samples
stages rather than tsinfer.infer(). This makes memory use and failure points easier to
inspect, and it allows mmap-backed ancestor generation on local Linux storage.

The input should be produced by ``f/arg.sh prep_vcf``. It must contain only
biallelic variants, phased GTs, and a normalized INFO/AA field. INFO/AA, rather
than VCF REF, is used as the ancestral state.

Large VCF-Zarr intermediates are created below a Linux temporary directory
(default: /tmp), not beside the final output on /mnt/d. This is important in
WSL because Zarr contains many small files and is extremely slow to create or
remove on a Windows-mounted filesystem.
"""



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


def infer_die(message: str) -> "NoReturn":
    raise SystemExit(message)


def infer_command_path(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        infer_die(f"Required command not found: {name}")
    return path


def infer_version_of(module: Any) -> str:
    return str(getattr(module, "__version__", "unknown"))


def infer_extract_sample_ids(vcf_path: Path) -> list[str]:
    """Read sample IDs with bcftools, avoiding Python gzip/header edge cases."""
    result = subprocess.run(
        [infer_command_path("bcftools"), "query", "-l", str(vcf_path)],
        check=True,
        text=True,
        capture_output=True,
    )
    samples = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if not samples:
        infer_die(f"No samples found in VCF: {vcf_path}")
    if len(samples) != len(set(samples)):
        infer_die(f"Duplicate sample IDs found in VCF: {vcf_path}")
    return samples


def infer_check_phasing(vcf_path: Path, max_records: int = 2000) -> dict[str, int]:
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
        infer_die(f"VCF has no records: {vcf_path}")
    if unphased_heterozygous:
        infer_die(
            f"VCF is not fully phased: found {unphased_heterozygous} unphased "
            f"heterozygous calls in the first {records} records"
        )
    return {
        "records_checked": records,
        "called_genotypes_checked": called,
        "heterozygous_genotypes_checked": heterozygous,
        "unphased_heterozygous_checked": unphased_heterozygous,
    }


def infer_scalar_string(value: Any) -> str:
    """Convert Zarr string/byte/scalar/length-one arrays to a plain string."""
    if hasattr(value, "tolist"):
        value = value.tolist()
    while isinstance(value, (list, tuple)) and len(value) == 1:
        value = value[0]
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value)


def infer_normalize_aa(value: Any) -> str:
    aa = infer_scalar_string(value).strip().upper()
    aa = aa.split("|", 1)[0].split(",", 1)[0]
    return aa


def infer_build_ancestral_array(zarr_path: Path, chunk_size: int = 250_000) -> tuple[str, int]:
    """Create a validated one-base ancestral-state array from variant_AA."""
    import numpy as np
    import zarr

    root = zarr.open(str(zarr_path), mode="r+")
    required = ("variant_AA", "variant_allele")
    missing = [name for name in required if name not in root]
    if missing:
        infer_die(
            "VCF-Zarr is missing "
            + ", ".join(missing)
            + ". The VCF must retain INFO/AA; REF is not an acceptable fallback."
        )
    aa_source = root["variant_AA"]
    alleles_source = root["variant_allele"]
    n_variants = int(aa_source.shape[0])
    if n_variants == 0:
        infer_die("VCF-Zarr contains no variants")
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
            aa = infer_normalize_aa(raw_aa)
            alleles = {
                infer_scalar_string(x).strip().upper()
                for x in allele_block[offset]
                if infer_scalar_string(x).strip() not in ("", ".")
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
        infer_die(
            "Invalid ancestral states remain after VCF preparation. Examples: "
            + "; ".join(invalid)
        )
    return out_name, n_variants


def infer_vcf_zarr_ploidy(zarr_path: Path) -> int:
    """Read the fixed GT ploidy dimension written by vcf2zarr."""
    import zarr

    root = zarr.open(str(zarr_path), mode="r")
    if "call_genotype" not in root:
        infer_die("VCF-Zarr is missing call_genotype")
    shape = root["call_genotype"].shape
    if len(shape) != 3 or int(shape[2]) < 1:
        infer_die(f"Unexpected call_genotype shape in VCF-Zarr: {shape}")
    return int(shape[2])


def infer_metadata_dict(obj: Any) -> dict[str, Any]:
    metadata = getattr(obj, "metadata", None)
    if isinstance(metadata, dict):
        return metadata
    if isinstance(metadata, bytes):
        try:
            return json.loads(metadata.decode())
        except Exception:
            return {}
    return {}


def infer_write_sample_map(tree_sequence: Any, sample_ids: list[str], output_path: Path) -> None:
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
                metadata = infer_metadata_dict(obj)
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


def infer_dump_atomic(tree_sequence: Any, output: Path) -> None:
    temp_output = output.with_name(f"{output.name}.part.{os.getpid()}")
    try:
        tree_sequence.dump(str(temp_output))
        os.replace(temp_output, output)
    finally:
        temp_output.unlink(missing_ok=True)


def infer_parse_args() -> argparse.Namespace:
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


def infer_main() -> None:
    args = infer_parse_args()
    if args.threads < 1:
        infer_die("--threads must be positive")
    if args.sample_match_min_work < 0:
        infer_die("--sample-match-min-work must be nonnegative")
    for name, value in (
        ("recombination rate", args.recombination_rate),
        ("mutation rate", args.mutation_rate),
        ("population size", args.population_size),
    ):
        if value <= 0:
            infer_die(f"{name} must be positive")
    args.vcf = args.vcf.resolve()
    if not args.vcf.is_file() or args.vcf.stat().st_size == 0:
        infer_die(f"Input VCF is missing or empty: {args.vcf}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    tmp_filesystem = subprocess.check_output(
        [infer_command_path("stat"), "-f", "-c", "%T", str(args.tmp_root)], text=True
    ).strip()
    if tmp_filesystem in {"v9fs", "9p", "drvfs"}:
        infer_die(
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

    sample_ids = infer_extract_sample_ids(args.vcf)
    phase_qc = infer_check_phasing(args.vcf, args.phase_check_records)
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
                [infer_command_path("vcf2zarr"), "convert", str(args.vcf), str(vcz_path)],
                check=True,
            )
        input_ploidy = infer_vcf_zarr_ploidy(vcz_path)
        x_identities = None
        if args.chr.upper().removeprefix('CHR') in ('X', '23'):

            haploid_path = work / f"{args.out.stem}.haploid.vcz"
            x_identities = x_haploid_zarr(vcz_path, haploid_path)
            vcz_path = haploid_path
            print(f"[ARG] Non-PAR X: {len(sample_ids)} individuals, {len(x_identities)} real haplotypes", flush=True)
        ancestral_array, n_input_variants = infer_build_ancestral_array(vcz_path)
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
            infer_dump_atomic(ancestors_ts, ancestors_ts_path)

        if inferred_path.is_file():
            print(f"[ARG] Reuse matched-samples checkpoint: {inferred_path}", flush=True)
            inferred_ts = tskit.load(str(inferred_path))
        elif args.sample_match_min_work == 0:
            print("[ARG] tsinfer.match_samples", flush=True)
            inferred_ts = tsinfer.match_samples(data, ancestors_ts, **match_kwargs)
            infer_dump_atomic(inferred_ts, inferred_path)
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
            infer_dump_atomic(inferred_ts, inferred_path)

        if x_identities is not None:

            inferred_ts = x_restore_individuals(inferred_ts, x_identities)
        if hasattr(tsdate, "preprocess_ts"):
            print("[ARG] tsdate.preprocess_ts", flush=True)
            processed_ts = tsdate.preprocess_ts(inferred_ts)
        else:
            processed_ts = inferred_ts
        infer_dump_atomic(processed_ts, processed_path)

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
            infer_die(
                "ARG sample-node count does not match VCF samples × ploidy "
                f"({dated_ts.num_samples} vs {len(sample_ids)} × {input_ploidy} = {expected_nodes})"
            )
        infer_dump_atomic(dated_ts, args.out)

        map_part = sample_map.with_name(f"{sample_map.name}.part.{os.getpid()}")
        infer_write_sample_map(dated_ts, sample_ids, map_part)
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
                "tsinfer": infer_version_of(tsinfer),
                "tsdate": infer_version_of(tsdate),
                "tskit": infer_version_of(tskit),
                "zarr": infer_version_of(zarr),
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





# ---- check ----
"""Fast structural validation for a completed tsinfer ARG and its sidecars."""



import argparse
import csv
import json
from pathlib import Path


from comm import enable_wide_csv_fields


enable_wide_csv_fields()


def check_fail(message: str) -> None:
    import tskit
    raise SystemExit(f"ERROR: ARG post-build check failed: {message}")


def check_resolved(value: str | Path) -> str:
    import tskit
    return str(Path(value).resolve())


def check_main() -> None:
    import tskit
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--qc-json", required=True, type=Path)
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--chr", required=True)
    args = parser.parse_args()

    try:
        ts = tskit.load(str(args.trees))
    except Exception as exc:
        check_fail(f"unreadable tree sequence {args.trees}: {exc}")
    if ts.num_trees < 1 or ts.num_sites < 1 or ts.num_samples < 1:
        check_fail(
            f"empty tree sequence: trees={ts.num_trees} sites={ts.num_sites} "
            f"samples={ts.num_samples}"
        )
    if ts.sequence_length <= 0:
        check_fail(f"invalid sequence length: {ts.sequence_length}")

    try:
        qc = json.loads(args.qc_json.read_text(encoding="utf-8"))
    except Exception as exc:
        check_fail(f"unreadable QC JSON {args.qc_json}: {exc}")
    expected_paths = {
        "input_vcf": check_resolved(args.input_vcf),
        "output_trees": check_resolved(args.trees),
        "sample_map": check_resolved(args.sample_map),
    }
    for key, expected in expected_paths.items():
        if check_resolved(qc.get(key, "")) != expected:
            check_fail(f"QC {key} does not match {expected}")
    expected_counts = {
        "sample_nodes": ts.num_samples,
        "sites": ts.num_sites,
        "mutations": ts.num_mutations,
        "nodes": ts.num_nodes,
        "edges": ts.num_edges,
        "trees": ts.num_trees,
    }
    for key, expected in expected_counts.items():
        if qc.get(key) != expected:
            check_fail(f"QC {key}={qc.get(key)!r}, tree sequence has {expected}")
    if float(qc.get("sequence_length", -1)) != float(ts.sequence_length):
        check_fail("QC sequence_length does not match tree sequence")
    if qc.get("ancestral_state_source") != "INFO/AA (VCF-Zarr variant_AA)":
        check_fail("QC does not record INFO/AA as the ancestral-state source")

    try:
        with args.sample_map.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["tree_node_id", "sample", "haplotype"]:
                check_fail(f"unexpected sample-map header: {reader.fieldnames}")
            rows = list(reader)
    except Exception as exc:
        check_fail(f"unreadable sample map {args.sample_map}: {exc}")
    if len(rows) != ts.num_samples:
        check_fail(f"sample-map rows={len(rows)}, expected {ts.num_samples}")
    try:
        mapped_nodes = [int(row["tree_node_id"]) for row in rows]
        if any(not row["sample"] or int(row["haplotype"]) < 1 for row in rows):
            check_fail("sample map contains an empty sample or invalid haplotype")
    except (TypeError, ValueError) as exc:
        check_fail(f"invalid sample-map value: {exc}")
    if len(mapped_nodes) != len(set(mapped_nodes)):
        check_fail("sample map contains duplicate tree node IDs")
    if set(mapped_nodes) != set(map(int, ts.samples())):
        check_fail("sample-map node IDs do not equal tree-sequence sample nodes")
    if args.chr == 'X' and 'sample_ploidies' not in qc:
        check_fail('X QC lacks biological sample ploidies; rebuild with the updated X workflow')
    if 'sample_ploidies' in qc:
        observed = {}
        for row in rows:
            observed.setdefault(row['sample'], []).append(int(row['haplotype']))
        expected = {name: list(range(1, int(ploidy) + 1)) for name, ploidy in qc['sample_ploidies'].items()}
        if {name: sorted(haps) for name, haps in observed.items()} != expected:
            check_fail('Sample map differs from biological sample ploidies')

    print(
        f"[ARG CHECK] PASS chr={args.chr} trees={ts.num_trees} "
        f"sites={ts.num_sites} samples={ts.num_samples} sequence_length={ts.sequence_length:g}"
    )





# ---- trace ----
"""Export a native ARG tree and sample map to GU TRACE's input contract."""


import argparse
import csv
import json
import os
from pathlib import Path



def trace_atomic_text(path: Path, text: str) -> None:
    import tskit
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, path)


def trace_read_map(path: Path, method: str) -> list[tuple[int, str, int]]:
    import tskit
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        rows = []
        if method == "needle":
            required = {"eid", "hap", "sample_node"}
            if not required <= fields:
                raise SystemExit(
                    f"Needle sample map requires {sorted(required)}: {path}"
                )
            for row in reader:
                rows.append(
                    (int(row["sample_node"]), str(row["eid"]), int(row["hap"]) + 1)
                )
        else:
            required = {"tree_node_id", "sample", "haplotype"}
            if not required <= fields:
                raise SystemExit(
                    f"tsinfer sample map requires {sorted(required)}: {path}"
                )
            for row in reader:
                rows.append(
                    (
                        int(row["tree_node_id"]),
                        str(row["sample"]),
                        int(row["haplotype"]),
                    )
                )
    return rows


def trace_validate(ts: tskit.TreeSequence, rows: list[tuple[int, str, int]], chrom: str) -> dict:
    import tskit
    if min(ts.num_trees, ts.num_samples, ts.num_sites) < 1:
        raise SystemExit(f"Empty tree sequence for chr{chrom}")
    nodes = list(map(int, ts.samples()))
    mapped = [row[0] for row in rows]
    if len(rows) != len(nodes) or set(mapped) != set(nodes) or len(set(mapped)) != len(mapped):
        raise SystemExit(
            f"TRACE sample map does not bijectively cover chr{chrom} sample nodes: "
            f"rows={len(rows)} samples={len(nodes)}"
        )
    positions = [site.position for site in ts.sites()]
    span = (max(positions) - min(positions)) / max(1.0, float(ts.sequence_length))
    if any(hap < 1 for _, _, hap in rows):
        raise SystemExit(f"TRACE haplotype numbers must be one-based for chr{chrom}")
    return {
        "chromosome": str(chrom),
        "trees": ts.num_trees,
        "nodes": ts.num_nodes,
        "samples": ts.num_samples,
        "sites": ts.num_sites,
        "mutations": ts.num_mutations,
        "sequence_length": ts.sequence_length,
        "site_span_fraction": span,
    }


def trace_write_tree(source: Path, dest: Path, ts: tskit.TreeSequence, method: str) -> str:
    import tskit
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f".{dest.name}.tmp.{os.getpid()}")
    tmp.unlink(missing_ok=True)
    offset = float(ts.metadata.get('offset', 0)) if method == 'needle' and isinstance(ts.metadata, dict) else 0
    if offset:
        # Native ARG-Needle coordinates are relative to the first input SNP.
        # GU loci and TRACE maps are chromosome coordinates, so translate once.
        tables = ts.dump_tables()
        tables.sequence_length += offset
        tables.edges.left += offset
        tables.edges.right += offset
        tables.sites.position += offset
        if len(tables.migrations):
            tables.migrations.left += offset
            tables.migrations.right += offset
        tables.metadata_schema = tskit.MetadataSchema.permissive_json()
        tables.metadata = {**ts.metadata, 'offset': 0, 'source_arg_offset': offset,
                           'coordinate_system': 'chromosome_bp'}
        tables.time_units = 'generations'
        tables.tree_sequence().dump(tmp)
        os.replace(tmp, dest)
        return 'chromosome_coordinate_copy'
    if ts.time_units == "generations":
        os.symlink(os.path.relpath(source.resolve(), dest.parent.resolve()), tmp)
        os.replace(tmp, dest)
        return "relative_symlink"
    if method != "needle" or ts.time_units not in ("unknown", ""):
        raise SystemExit(
            f"TRACE requires node times in generations; got {ts.time_units!r}: {source}"
        )
    tables = ts.dump_tables()
    tables.time_units = "generations"
    tables.provenances.add_row(
        record=json.dumps(
            {
                "schema_version": "1.0.0",
                "software": {"name": "arg.trace", "version": "1"},
                "parameters": {
                    "source_method": "needle",
                    "operation": "declare_ARG_Needle_normalized_node_times_as_generations",
                },
            },
            separators=(",", ":"),
        )
    )
    tables.tree_sequence().dump(tmp)
    os.replace(tmp, dest)
    return "converted_copy"


def trace_main() -> None:
    import tskit
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, choices=("needle", "tsinfer", "threads"))
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--chr", required=True)
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    for path in (args.tree, args.sample_map):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing ARG input: {path}")
    out_tree = args.out_dir / f"chr{args.chr}.trees"
    out_map = args.out_dir / f"chr{args.chr}.sample_map.tsv"
    out_qc = args.out_dir / f"chr{args.chr}.arg_qc.json"
    out_meta = args.out_dir / f"chr{args.chr}.trace.meta.json"
    source_signature = {
        "tree": [str(args.tree.resolve()), args.tree.stat().st_size, args.tree.stat().st_mtime_ns],
        "sample_map": [
            str(args.sample_map.resolve()),
            args.sample_map.stat().st_size,
            args.sample_map.stat().st_mtime_ns,
        ],
        "method": args.method,
        "format": "trace",
    }
    outputs = (out_tree, out_map, out_qc, out_meta)
    if args.check:
        for path in outputs:
            if not path.exists() or path.stat().st_size == 0:
                raise SystemExit(f"Missing TRACE-format ARG output: {path}")
        old = json.loads(out_meta.read_text())
        if old.get("source") != source_signature:
            raise SystemExit(
                f"TRACE-format ARG source provenance differs for chr{args.chr}; "
                "rerun build with --format trace"
            )
        out_ts = tskit.load(out_tree)
        rows = trace_read_map(out_map, "tsinfer")
        stats = trace_validate(out_ts, rows, args.chr)
        if out_ts.time_units != "generations":
            raise SystemExit(f"TRACE tree time units are not generations: {out_tree}")
        print(
            f"TRACE format CHECK PASS method={args.method} chr{args.chr} "
            f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
        )
        return

    if not args.replace and all(path.exists() and path.stat().st_size > 0 for path in outputs):
        try:
            old = json.loads(out_meta.read_text())
            if old.get("source") == source_signature:
                out_ts = tskit.load(out_tree)
                rows = trace_read_map(out_map, "tsinfer")
                stats = trace_validate(out_ts, rows, args.chr)
                if out_ts.time_units == "generations":
                    print(
                        f"TRACE format SKIP method={args.method} chr{args.chr} "
                        f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
                    )
                    return
        except (OSError, ValueError, json.JSONDecodeError, tskit.FileFormatError):
            pass

    ts = tskit.load(args.tree)
    rows = trace_read_map(args.sample_map, args.method)
    stats = trace_validate(ts, rows, args.chr)
    storage = trace_write_tree(args.tree, out_tree, ts, args.method)
    map_text = "tree_node_id\tsample\thaplotype\n" + "".join(
        f"{node}\t{sample}\t{hap}\n" for node, sample, hap in rows
    )
    trace_atomic_text(out_map, map_text)
    out_ts = tskit.load(out_tree)
    stats = trace_validate(out_ts, trace_read_map(out_map, "tsinfer"), args.chr)
    if out_ts.time_units != "generations":
        raise SystemExit(f"TRACE tree time units are not generations: {out_tree}")
    qc = {
        **stats,
        "method": args.method,
        "format": "trace",
        "time_units": out_ts.time_units,
        "storage": storage,
    }
    trace_atomic_text(out_qc, json.dumps(qc, indent=2, sort_keys=True) + "\n")
    trace_atomic_text(
        out_meta,
        json.dumps({"source": source_signature, "output": qc}, indent=2, sort_keys=True)
        + "\n",
    )
    print(
        f"TRACE format PASS method={args.method} chr{args.chr} storage={storage} "
        f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
    )





# ---- unique ----
"""Drop all records at repeated positions from sorted PLINK Oxford HAPS.

Keep only one pending row in memory. Do not select an arbitrary allele from
split multiallelic sites. Replace the input only after a successful full scan.
"""
import argparse
import json
import os
from pathlib import Path


def unique_clean(path, chromosome):
    path = Path(path)
    temporary = path.with_name(path.name + '.unique.next')
    total = kept = duplicate_sites = duplicate_records = 0
    previous = None
    pending = None
    count = 0
    try:
        with path.open() as source, temporary.open('w', newline='\n') as output:
            for number, line in enumerate(source, 1):
                if not line.strip():
                    continue
                fields = line.split(maxsplit=5)
                if len(fields) != 6 or fields[0].removeprefix('chr') != str(chromosome):
                    raise ValueError(f'HAPS line {number}: malformed row or wrong chromosome')
                position = int(fields[2])
                if position < 1 or (previous is not None and position < previous):
                    raise ValueError(f'HAPS line {number}: invalid/unsorted position {position} after {previous}')
                total += 1
                if position == previous:
                    if count == 1:
                        duplicate_sites += 1
                        duplicate_records += 1
                    duplicate_records += 1
                    count += 1
                    pending = None
                else:
                    if pending is not None:
                        output.write(pending)
                        kept += 1
                    previous, pending, count = position, line, 1
            if pending is not None:
                output.write(pending)
                kept += 1
        if kept < 2:
            raise ValueError(f'Too few unique-position HAPS variants: {kept}')
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return dict(input_variants=total, retained_variants=kept,
                duplicate_positions=duplicate_sites, removed_variants=duplicate_records,
                policy='exclude-all-repeated-positions')


def unique_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--haps', required=True)
    parser.add_argument('--chr', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    result = unique_clean(args.haps, args.chr)
    Path(args.out).write_text(json.dumps(result, indent=2) + '\n')
    print('HAPS position QC: ' + json.dumps(result))






def cache_main():
    #!/usr/bin/env python3
    """Stable request fingerprints. Files use canonical names and nanosecond stat data."""
    import argparse
    import hashlib
    import json
    from pathlib import Path

    p = argparse.ArgumentParser()
    p.add_argument('--file', action='append', default=[])
    p.add_argument('--text-file', action='append', default=[])
    p.add_argument('--value', action='append', default=[])
    a = p.parse_args()
    files = []
    for name in a.file + a.text_file:
        path = Path(name).resolve(strict=True)
        s = path.stat()
        row = [str(path), s.st_size, s.st_mtime_ns]
        if name in a.text_file:
            row.append(hashlib.sha256(path.read_bytes()).hexdigest())
        files.append(row)
    print(json.dumps({'schema': 2, 'values': a.value, 'files': files}, sort_keys=True))


def validate_vcf_main():
    #!/usr/bin/env python3
    """Full, streaming phase/AA/position validation before committing prepared VCF."""
    import argparse
    import pysam

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vcf')
    parser.add_argument('--chr', default='')
    parser.add_argument('--build', choices=('37', '38'), default='37')
    parser.add_argument('--x-psam')
    args = parser.parse_args()
    with pysam.VariantFile(args.vcf) as vcf:
        ploidies = None
        if args.chr == 'X':
            if args.x_psam:
                sexes = x_read_sexes(args.x_psam)
                if not set(vcf.header.samples) <= sexes.keys():
                    raise SystemExit('X VCF contains samples absent from PSAM')
                ploidies = {n: 1 if sexes[n] == 1 else 2 for n in vcf.header.samples}
        previous = None
        count = 0
        for r in vcf:
            count += 1
            if previous and (r.contig != previous[0] or r.pos <= previous[1]):
                raise SystemExit(f'Non-unique/unsorted positions: {r.contig}:{r.pos}')
            previous = (r.contig, r.pos)
            if args.chr == 'X':
                if r.contig.removeprefix('chr') not in ('X','23') or not x_nonpar(r.pos, args.build):
                    raise SystemExit(f'Expected non-PAR X: {r.contig}:{r.pos}')
                if ploidies is None:
                    ploidies = {n: len(c['GT']) for n, c in r.samples.items()}
            aa = r.info.get('AA')
            if len(r.alleles) != 2 or aa not in r.alleles or any(x not in ('A','C','G','T') for x in r.alleles):
                raise SystemExit(f'Invalid alleles/AA: {r.contig}:{r.pos}')
            for name, call in r.samples.items():
                gt = call.get('GT')
                if gt is None or len(gt) not in (1,2): raise SystemExit(f'Invalid ploidy at {r.pos}: {name}')
                if ploidies is not None and len(gt) != ploidies[name]:
                    raise SystemExit(f'Incorrect/changing non-PAR X ploidy at {r.pos}: {name}; expected {ploidies[name]}, got {len(gt)}. Supply a VCF with haploid males and diploid females.')
                if any(x is not None and x not in (0,1) for x in gt): raise SystemExit('Invalid GT allele')
                if len(gt)==2 and None not in gt and gt[0]!=gt[1] and not call.phased:
                    raise SystemExit(f'Unphased heterozygote at {r.contig}:{r.pos}: {name}')
        if count == 0: raise SystemExit('No variants')
    print(f'Validated all {count} VCF records: phased GT, valid AA, unique positions')


COMMANDS = {'cache': cache_main, 'validate-vcf': validate_vcf_main,
            'unique-haps': unique_main, 'needle-x': needle_x_main,
            'infer-tsinfer': infer_main, 'check-tsinfer': check_main, 'trace': trace_main}

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print('Usage: arg.py {' + ','.join(COMMANDS) + '} [options]')
        return
    command = sys.argv.pop(1)
    if command not in COMMANDS:
        raise SystemExit('Unknown ARG command: ' + command)
    COMMANDS[command]()

if __name__ == '__main__':
    main()
