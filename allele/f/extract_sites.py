#!/usr/bin/env python3
"""Extract observed SNV loci and expand each locus to its three possible ALT bases.

The denominator is conditional on loci observed as variable in the input VCF(s).
For each locus, observed alternatives come from the VCF; absent alternatives are
counterfactual candidates for AlphaGenome. Coordinates are GRCh38 1-based VCF
coordinates. REF is checked against FASTA when possible.
"""
from __future__ import annotations
import argparse
import json
import os
import re
from collections import Counter
from contextlib import ExitStack
from heapq import merge
from itertools import groupby
from pathlib import Path
from parquet_io import ParquetSink
from vcf_inputs import chrom_key, resolve_inputs, validate_inputs
import numpy as np
import pandas as pd
import pysam
from cyvcf2 import VCF
from tqdm import tqdm

BASES = np.array(list("ACGT"))


def parse_args():
    p = argparse.ArgumentParser()
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--vcf-dir", help="directory of VCFs from one cohort")
    source.add_argument("--vcf", help="one VCF file, including local cohort VCFs")
    p.add_argument(
        "--region", help="CHR:START-END, 1-based inclusive; requires indexed --vcf"
    )
    p.add_argument("--chr", help="autosomes, all, or 1,2,... (directory input)")
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--batch-sites", type=int, default=10000)
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--sites-out", required=True)
    p.add_argument(
        "--include-filtered", action="store_true", help="include non-PASS records"
    )
    args = p.parse_args()
    if args.batch_sites < 1:
        p.error("--batch-sites must be positive")
    if args.region:
        if not args.vcf:
            p.error("--region requires --vcf (select one indexed file)")
        match = re.fullmatch(r"([^:\s,]+):([1-9][0-9]*)-([1-9][0-9]*)", args.region)
        if not match or int(match[2]) > int(match[3]):
            p.error("--region must be CHR:START-END with 1 <= START <= END")
        args.region_bounds = (match[1], int(match[2]), int(match[3]))
    return args


def info_array(v, key, n):
    try:
        val = v.INFO.get(key)
    except Exception:
        val = None
    if val is None:
        return [np.nan] * n
    if isinstance(val, (tuple, list, np.ndarray)):
        arr = list(val)
    else:
        arr = [val]
    arr = (arr + [np.nan] * n)[:n]
    return [float(x) if x is not None else np.nan for x in arr]


def normalize_chrom(c, fasta_refs):
    if c in fasta_refs:
        return c
    c2 = c.removeprefix("chr")
    if "chr" + c2 in fasta_refs:
        return "chr" + c2
    if c2 in fasta_refs:
        return c2
    return c


def locus_rows(key, d, fa, qc):
    chrom, pos, ref = key
    cand, sites = [], []
    obs = d["alts"]
    occ = len(obs)
    if occ < 1 or occ > 3:
        return [], None
    try:
        context = fa.fetch(chrom, max(0, pos - 2), pos + 1).upper()
    except Exception:
        context = "NNN"
    if len(context) != 3:
        context = "NNN"
    afs = np.array([x["af"] for x in obs.values() if not np.isnan(x["af"])], float)
    if afs.sum() > 1 + 1e-6:
        raise SystemExit(
            f"ALT frequencies sum to >1 at {chrom}:{pos}; check overlapping inputs/cohorts"
        )
    ref_af = max(0.0, 1.0 - float(afs.sum())) if len(afs) == occ else np.nan
    if np.isfinite(ref_af) and ref_af <= 1e-8:
        qc["reference_absent_loci"] += 1
        return [], None
    probs = [ref_af] + [x["af"] for x in obs.values() if not np.isnan(x["af"])]
    probs = np.array([p for p in probs if p > 0 and np.isfinite(p)], float)
    entropy = (
        float(-(probs * np.log(probs)).sum() / np.log(max(2, len(probs))))
        if len(probs) > 1
        else 0.0
    )
    if not np.isfinite(ref_af):
        entropy = np.nan
    site_id = f"{chrom}:{pos}:{ref}"
    primary = max(
        sorted(obs),
        key=lambda x: (-np.inf if np.isnan(obs[x]["af"]) else obs[x]["af"]),
    )
    sites.append(
        dict(
            site_id=site_id,
            chrom=chrom,
            pos=pos,
            ref=ref,
            context=context,
            occupancy=occ,
            allele_count=occ + 1,
            entropy=entropy,
            primary_alt=primary,
            ref_af=ref_af,
        )
    )
    for alt in BASES[BASES != ref]:
        alt = str(alt)
        o = alt in obs
        m = obs.get(alt, {"af": 0.0, "ac": 0.0, "an": np.nan})
        is_transition = (ref + alt) in {"AG", "GA", "CT", "TC"}
        cpg = (
            len(context) == 3 and context[1] == "C" and context[2] == "G" and alt == "T"
        ) or (
            len(context) == 3
            and context[0] == "C"
            and context[1] == "G"
            and ref == "G"
            and alt == "A"
        )
        cand.append(
            dict(
                site_id=site_id,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                variant_id=f"{chrom}:{pos}:{ref}>{alt}",
                observed=int(o),
                occupancy=occ,
                allele_count=occ + 1,
                af=float(m["af"]) if np.isfinite(m["af"]) else np.nan,
                ac=float(m["ac"]) if np.isfinite(m["ac"]) else np.nan,
                an=float(m["an"]) if np.isfinite(m["an"]) else np.nan,
                primary_alt=primary,
                is_primary=int(alt == primary),
                context=context,
                ref_alt=ref + ">" + alt,
                transition=int(is_transition),
                cpg_transition=int(cpg),
                entropy=entropy,
                ref_af=ref_af,
            )
        )
    return cand, sites[0]


def read_records(path, reader, args, fasta_refs, qc):
    if args.region:
        region_chrom, start, end = args.region_bounds
        region_chrom = normalize_chrom(region_chrom, set(reader.seqnames))
        if region_chrom not in reader.seqnames:
            raise ValueError(f"Region chromosome is absent from VCF: {args.region}")
        records = reader(f"{region_chrom}:{start}-{end}")
    else:
        records = reader
    previous = None
    for v in tqdm(records, desc=path.name, unit="var", mininterval=5):
        if args.region and not start <= v.POS <= end:
            continue
        qc["records"] += 1
        chrom = normalize_chrom(v.CHROM, fasta_refs)
        order = (chrom_key(chrom), int(v.POS))
        if previous is not None and order < previous:
            raise ValueError(
                f"VCF must be coordinate sorted: {path} at {chrom}:{v.POS}"
            )
        previous = order
        if len(v.REF) != 1 or v.REF not in "ACGT":
            continue
        alts = [
            (i, x)
            for i, x in enumerate(v.ALT or [])
            if x and len(x) == 1 and x in "ACGT" and x != v.REF
        ]
        if not alts:
            continue
        if not args.include_filtered and v.FILTER not in (None, "PASS", "."):
            qc["filtered_records"] += 1
            continue
        af = info_array(v, "AF", len(v.ALT))
        ac = info_array(v, "AC", len(v.ALT))
        try:
            an = float(v.INFO.get("AN"))
        except (TypeError, ValueError):
            an = np.nan
        # Copy INFO values: cyvcf2 record storage must not escape the iterator.
        yield order, chrom, int(v.POS), v.REF, alts, af, ac, an


def main():
    a = parse_args()
    try:
        files = resolve_inputs(a.vcf, a.vcf_dir, a.chr)
        validate_inputs(files, require_index=bool(a.region))
    except ValueError as error:
        raise SystemExit(str(error))
    if not os.path.isfile(a.fasta + ".fai"):
        raise SystemExit(f"Missing FASTA index; run: samtools faidx {a.fasta}")
    with pysam.FastaFile(a.fasta) as fa:
        if a.check_only:
            print(
                f"Inputs ready: {len(files)} VCF files; FASTA index opens successfully."
            )
            print(
                "Checks: requested files/indexes, pending .aria2, official gnomAD sizes, BGZF end markers."
            )
            print("Full-file MD5 is checked by gnomad_download.sh, not repeated here.")
            for path in files:
                print(path)
            return
        for path in [a.out, a.sites_out]:
            if os.path.exists(path):
                raise SystemExit(
                    f"Output already exists: {path}; use a new --out directory to preserve matching downstream results"
                )
        qc = Counter()
        inputs = []
        counts = Counter()
        sinks = [ParquetSink(a.out), ParquetSink(a.sites_out)]
        candidates, sites = [], []
        try:
            with ExitStack() as stack:
                streams = []
                for path in files:
                    reader = VCF(str(path), lazy=True)
                    stack.callback(reader.close)
                    inputs.append(
                        dict(
                            path=str(path.resolve()),
                            bytes=path.stat().st_size,
                            sample_count=len(reader.samples),
                        )
                    )
                    reader.set_samples([])
                    streams.append(
                        read_records(path, reader, a, set(fa.references), qc)
                    )
                records = merge(*streams, key=lambda r: r[0])
                for _, group in groupby(records, key=lambda r: r[0]):
                    d = {"alts": {}}
                    key = None
                    for _, chrom, pos, ref, alts, af, ac, an in group:
                        if key and key[2] != ref:
                            raise ValueError(f"Conflicting REF bases at {chrom}:{pos}")
                        key = (chrom, pos, ref)
                        for i, alt in alts:
                            old = d["alts"].get(alt, {"af": np.nan})
                            af_i, ac_i = af[i], ac[i]
                            if (
                                np.isnan(af_i)
                                and np.isfinite(ac_i)
                                and np.isfinite(an)
                                and an > 0
                            ):
                                af_i = ac_i / an
                            if (np.isfinite(ac_i) and ac_i <= 0) or (
                                np.isfinite(af_i) and af_i <= 0
                            ):
                                qc["zero_count_alternatives"] += 1
                                continue
                            if np.isfinite(af_i) and not 0 <= af_i <= 1:
                                raise ValueError(f"Invalid AF at {chrom}:{pos}: {af_i}")
                            if np.isnan(old["af"]) or (
                                np.isfinite(af_i) and af_i > old["af"]
                            ):
                                d["alts"][alt] = {"af": af_i, "ac": ac_i, "an": an}
                    chrom, pos, ref = key
                    try:
                        ref_fa = fa.fetch(chrom, pos - 1, pos).upper()
                    except (KeyError, ValueError):
                        ref_fa = None
                    if ref_fa != ref:
                        qc["reference_mismatch"] += 1
                        raise ValueError(
                            f"SNV records disagree with FASTA at {chrom}:{pos}:{ref}; check GRCh38 build/contig names"
                        )
                    new_candidates, site = locus_rows(key, d, fa, qc)
                    if site is None:
                        continue
                    counts[site["occupancy"]] += 1
                    candidates.extend(new_candidates)
                    sites.append(site)
                    if len(sites) >= a.batch_sites:
                        sinks[0].write(pd.DataFrame(candidates))
                        sinks[1].write(pd.DataFrame(sites))
                        candidates, sites = [], []
                if not counts:
                    raise ValueError(
                        "No eligible observed SNV loci remain after filtering; check VCF, counts and FASTA"
                    )
                if sites:
                    sinks[0].write(pd.DataFrame(candidates))
                    sinks[1].write(pd.DataFrame(sites))
                for sink in sinks:
                    sink.commit()
        except ValueError as error:
            raise SystemExit(str(error))
        finally:
            for sink in sinks:
                sink.close()
            Path(a.out + ".qc.json").write_text(json.dumps(dict(qc), indent=2) + "\n")
    with open(a.out + ".inputs.json", "w") as handle:
        json.dump(
            dict(
                files=inputs,
                region=a.region,
                chromosomes=a.chr,
                counts_source="VCF INFO AC/AN/AF; genotypes not recalculated",
                site_count=sum(counts.values()),
                occupancy_counts=dict(counts),
            ),
            handle,
            indent=2,
        )
    print(
        f"loci={sum(counts.values()):,} candidates={3*sum(counts.values()):,} biallelic={counts[1]:,} triallelic={counts[2]:,} tetra={counts[3]:,}"
    )


if __name__ == "__main__":
    main()
