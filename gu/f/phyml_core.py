#!/usr/bin/env python3
"""Shared VCF and haplotype readers for the GU GWAS risk-core workflow."""
from __future__ import annotations

import csv
import re
import subprocess
import tempfile
from array import array
from collections import Counter
from typing import Iterator
from functools import lru_cache
from pathlib import Path

from comm import enable_wide_csv_fields


enable_wide_csv_fields()

BASES = {"A", "C", "G", "T"}


# 🚩 Locus errors and command execution
class SkipLocus(RuntimeError):
    """A chromosome-window has no usable data and may be skipped."""

    def __init__(self, message: str, code: str = "low_sequence_information", details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.details = details or {}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def run(args: list[str]) -> str:
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode:
        fail(f"command failed ({p.returncode}): {' '.join(args)}\n{p.stderr.strip()}")
    return p.stdout


# 🚩 Sample metadata and VCF lookup
def read_sexes(path: Path) -> dict[str, str]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample" not in reader.fieldnames or "sex" not in reader.fieldnames:
            fail(f"{path} must contain tab-separated sample and sex columns")
        out = {}
        for row in reader:
            sex = str(row["sex"]).strip().lower()
            if sex in {"1", "m", "male"}: sex = "male"
            elif sex in {"2", "f", "female"}: sex = "female"
            elif sex in {"", "0", "na", "n/a", ".", "unknown", "u"}: sex = "unknown"
            else: fail(f"invalid sex for sample {row['sample']}: {row['sex']}")
            out[row["sample"]] = sex
    return out


@lru_cache(maxsize=None)
def vcf_path(vcf_dir: Path, chrom: str) -> Path:
    candidates = [
        vcf_dir / f"chr{chrom}.vcf.gz",
        vcf_dir / f"{chrom}.vcf.gz",
        vcf_dir / f"chr{chrom}.bcf",
        vcf_dir / f"{chrom}.bcf",
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size:
            return path
    fail(f"modern VCF for chr{chrom} not found under {vcf_dir}")


@lru_cache(maxsize=None)
def contigs(vcf: Path) -> dict[str, int | None]:
    header = run(["bcftools", "view", "-h", str(vcf)])
    out: dict[str, int | None] = {}
    for line in header.splitlines():
        m = re.match(r"##contig=<ID=([^,>]+)(?:,length=([0-9]+))?", line, re.I)
        if m:
            out[m.group(1)] = int(m.group(2)) if m.group(2) else None
    return out


@lru_cache(maxsize=None)
def vcf_contig(vcf: Path, chrom: str) -> str:
    names = contigs(vcf)
    for candidate in (chrom, f"chr{chrom}", "23" if chrom == "X" else chrom):
        if candidate in names:
            return candidate
    fail(f"VCF {vcf} has no contig for chr{chrom}")


@lru_cache(maxsize=None)
def archaic_vcf(root: Path, ref: str, chrom: str) -> Path:
    aliases = [ref]
    if ref.lower() == "denisova":
        aliases += ["Denisovan"]
    dirs = []
    for alias in aliases:
        dirs += [root / alias, root / "avcf" / alias]
    patterns = [f"*chr{chrom}_*.vcf.gz", f"*chr{chrom}.*.vcf.gz", f"*chr{chrom}.vcf.gz", f"*{chrom}.vcf.gz"]
    for directory in dirs:
        if not directory.is_dir():
            continue
        for pattern in patterns:
            found = sorted(x for x in directory.glob(pattern) if x.is_file() and x.stat().st_size)
            if found:
                return found[0]
    fail(f"archaic VCF ref={ref} chr={chrom} not found under {root}")


@lru_cache(maxsize=None)
def vcf_samples(vcf: Path) -> tuple[str, ...]:
    return tuple(x for x in run(["bcftools", "query", "-l", str(vcf)]).splitlines() if x)


# 🚩 Streaming VCF queries
def query_rows(vcf: Path, region: str, include_aa: bool = False) -> tuple[list[str], Iterator[list[str]]]:
    samples = list(vcf_samples(vcf))
    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT" + ("\t%INFO/AA" if include_aa else "") + "[\t%GT]\n"
    command = ["bcftools", "query", "-r", region, "-f", fmt, str(vcf)]
    def stream():
        # Do not retain the full sample-by-variant text or millions of GT strings.
        # A file for stderr avoids pipe deadlocks when bcftools emits warnings.
        with tempfile.TemporaryFile(mode="w+t") as errors:
            proc = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=errors)
            try:
                for line in proc.stdout:
                    if line.strip():
                        yield line.rstrip("\r\n").split("\t")
                if proc.wait():
                    errors.seek(0)
                    fail(f"command failed ({proc.returncode}): {' '.join(command)}\n{errors.read().strip()}")
            finally:
                proc.stdout.close()
                if proc.poll() is None:
                    proc.terminate()
                proc.wait()
    return samples, stream()


# 🚩 Genotype and ancestral allele decoding
def called_base(ref: str, alt: str, gt: str, phased_required: bool) -> tuple[str, str]:
    alleles = [ref] + alt.split(",")
    if gt in {".", "./.", ".|."}:
        return "N", "N"
    if phased_required and "/" in gt and gt.split("/")[0] != gt.split("/")[-1]:
        return "N", "N"
    parts = re.split(r"[/|]", gt)
    if len(parts) == 1:
        parts *= 2
    if len(parts) != 2 or any(not x.isdigit() or int(x) >= len(alleles) for x in parts):
        return "N", "N"
    return alleles[int(parts[0])].upper(), alleles[int(parts[1])].upper()


def called_haploid_base(ref: str, alt: str, gt: str) -> str:
    """Return one male-X allele and reject every diploid representation."""
    alleles = [ref] + alt.split(",")
    if gt in {".", "", "./.", ".|."}:
        return "N"
    parts = re.split(r"[/|]", gt)
    if len(parts) != 1 or not parts[0].isdigit() or int(parts[0]) >= len(alleles):
        fail(f"chrX.male contains a non-haploid/heterozygous GT: {gt}")
    return alleles[int(parts[0])].upper()


def ancestral_allele(value: str, ref: str, alt: str, allow_third_allele: bool = False) -> str:
    """Normalize 1KG INFO/AA values such as A, A|, or A|||."""
    token = re.split(r"[|,]", str(value).strip().upper(), maxsplit=1)[0]
    return token if len(token) == 1 and token in BASES and (allow_third_allele or token in {ref, alt}) else "N"


def alt_dosage(gt: str, ploidy: int) -> int | None:
    """Return ALT-allele dosage without requiring phased genotypes."""
    if gt in {".", "", "./.", ".|."}:
        return None
    parts = re.split(r"[/|]", gt)
    if len(parts) != ploidy or any(not x.isdigit() or int(x) > 1 for x in parts):
        return None
    return sum(int(x) for x in parts)


# 🚩 Modern haplotypes
def modern_data(
    vcf: Path,
    contig: str,
    locus: dict,
    min_mac: int,
    sexes: dict[str, str],
    x_male_only: bool,
    x_par_diploid: bool,
    ancestral_any_base: bool = False,
):
    samples, rows = query_rows(vcf, f"{contig}:{locus['start'] + 1}-{locus['end']}", include_aa=True)
    if not samples:
        fail(f"modern VCF has no samples: {vcf}")
    missing = [sample for sample in samples if sample not in sexes] if sexes else []
    if sexes and missing:
        fail(f"samples absent from samples.txt (first examples): {','.join(missing[:10])}")
    sites = []
    anchors = []
    coordinate_anchor = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)", locus["name"], re.I)
    if locus["chrom"] == "X" and x_male_only and sexes:
        nonmale = [sample for sample in samples if sexes.get(sample) != "male"]
        if nonmale:
            fail(
                "chrX.male VCF contains samples not marked male "
                f"(first examples): {','.join(nonmale[:10])}"
            )
    for fields in rows:
        if len(fields) < 6 + len(samples):
            continue
        _, pos, vid, ref, alt, aa, *gts = fields
        pos_i = int(pos)
        coordinate_match = coordinate_anchor and (coordinate_anchor[1], int(coordinate_anchor[2]), coordinate_anchor[3].upper(), coordinate_anchor[4].upper()) == (locus["chrom"], pos_i, ref.upper(), alt.upper())
        if locus["name"] in vid.split(";") or coordinate_match:
            anchors.append((pos_i, ref, alt))
        if len(ref) != 1 or len(alt) != 1 or ref.upper() not in BASES or alt.upper() not in BASES:
            continue
        haps = []
        ploidy = 1 if locus["chrom"] == "X" and x_male_only else 2
        # GT vocabulary is tiny; decode once per site, not once per sample.
        decoded = {}
        dosage_values = array('b')
        for sample, gt in zip(samples, gts):
            if gt not in decoded:
                if locus["chrom"] == "X" and x_male_only:
                    pair = (called_haploid_base(ref.upper(), alt.upper(), gt), "N")
                elif locus["chrom"] != "X" or x_par_diploid:
                    pair = called_base(ref.upper(), alt.upper(), gt, phased_required=True)
                else:
                    fail("chrX input requires an explicit male non-PAR or diploid PAR mode")
                dose = alt_dosage(gt, ploidy)
                decoded[gt] = (pair, -1 if dose is None else dose)
            pair, dose = decoded[gt]
            haps.extend(pair)
            dosage_values.append(dose)
        haps = ''.join(haps)
        counts = Counter(x for x in haps if x in BASES)
        if len(counts) < 2 or min(counts.values()) < min_mac:
            continue
        sites.append(dict(pos=pos_i, vid=vid, ref=ref.upper(), alt=alt.upper(),
                          ancestral=ancestral_allele(aa,ref.upper(),alt.upper(),ancestral_any_base),haps=haps,dosages=dosage_values))
    if not sites:
        raise SkipLocus(f"no polymorphic biallelic SNPs with MAC >= {min_mac}")
    target_pos = anchors[0][0] if anchors else (locus["core_start"] + locus["core_end"]) // 2
    # LD selection needs an anchor that survived the biallelic/MAC filters.  An
    # exact named anchor is preferred; otherwise use the closest eligible SNP
    # and report its actual position downstream.
    anchor_pos = min(sites, key=lambda x: (abs(x["pos"] - target_pos), x["pos"]))["pos"]
    return samples, sites, anchor_pos


# 🚩 Archaic haplotypes
def archaic_calls(vcf: Path, contig: str, locus: dict, modern_sites: list[dict], allow_third_allele: bool = False) -> dict[int, str]:
    _, rows = query_rows(vcf, f"{contig}:{locus['start'] + 1}-{locus['end']}")
    wanted = {x["pos"]: x for x in modern_sites}
    calls: dict[int, str] = {}
    for fields in rows:
        if len(fields) < 6:
            continue
        _, pos, _, ref, alt, *gts = fields
        pos_i = int(pos)
        if pos_i not in wanted:
            continue
        observed = []
        for gt in gts:
            observed.extend(called_base(ref.upper(), alt.upper(), gt, phased_required=False))
        # A heterozygote containing one allele outside the modern pair is
        # still a heterozygote; filtering that allele first creates a false
        # homozygous match. Missing/heterozygous archaic calls remain unknown.
        calls[pos_i] = observed[0] if (observed and len(set(observed)) == 1
            and observed[0] in BASES and (allow_third_allele or observed[0] in {wanted[pos_i]["ref"], wanted[pos_i]["alt"]})) else "N"
    return calls


