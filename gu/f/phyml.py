#!/usr/bin/env python3
"""Locus-based modern/archaic haplotype comparison and PhyML workflow for GU."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path


BASES = {"A", "C", "G", "T"}

CHROM_LENGTHS = {
    "37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560,
    },
    "38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
        "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
        "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
        "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895,
    },
}

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


def norm_chr(value: str) -> str:
    value = re.sub(r"^chr", "", str(value), flags=re.I).upper()
    return "X" if value == "23" else value


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return value or "locus"


def read_locus_map(path: Path | None) -> dict[tuple[str, int, int, str], tuple[int, int]]:
    if path is None:
        return {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chr", "core_start", "core_end", "analysis_start", "analysis_end", "locus_id"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            fail(f"{path} is not a GU loci map")
        return {
            (
                norm_chr(row["chr"]),
                int(row["analysis_start"]),
                int(row["analysis_end"]),
                row["locus_id"],
            ): (int(row["core_start"]), int(row["core_end"]))
            for row in reader
        }


def read_loci(path: Path, map_path: Path | None = None) -> list[dict]:
    core_map = read_locus_map(map_path)
    rows = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                fail(f"{path}:{line_no}: expected at least 3 BED columns")
            chrom, start, end = norm_chr(fields[0]), int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 else f"locus_{chrom}_{start}_{end}"
            if start < 0 or end <= start:
                fail(f"{path}:{line_no}: invalid interval {start}-{end}")
            core_start, core_end = core_map.get((chrom, start, end, name), (start, end))
            rows.append(dict(chrom=chrom, start=start, end=end, name=name,
                             core_start=core_start, core_end=core_end))
    if not rows:
        fail(f"no loci in {path}")
    seen: Counter[str] = Counter()
    for row in rows:
        base = safe_name(row["name"])
        seen[base] += 1
        row["locus_id"] = base if seen[base] == 1 else f"{base}_{seen[base]}"
    return rows


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


def annotate_loci(loci: list[dict], x_par_diploid: bool = False) -> list[dict]:
    """Label X loci after GU has selected either chrX.male or chrXY."""
    return [
        {
            **locus,
            "x_class": (
                "par_diploid" if locus["chrom"] == "X" and x_par_diploid
                else "male_nonpar" if locus["chrom"] == "X"
                else "autosome"
            ),
        }
        for locus in loci
    ]


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


def chromosome_loci(chrs: list[str], vcf_dir: Path, window_bp: int, build: str) -> list[dict]:
    rows = []
    for chrom in chrs:
        vcf = vcf_path(vcf_dir, chrom)
        contig = vcf_contig(vcf, chrom)
        length = contigs(vcf).get(contig)
        if not length:
            length = CHROM_LENGTHS[build].get(chrom)
            print(f"PHYML CHECK: chr{chrom} VCF header has no contig length; using GRCh{build} length={length}", flush=True)
        if not length:
            fail(f"--chr phyml needs a contig length in the VCF header: {vcf} ({contig})")
        for start in range(0, length, window_bp):
            end = min(length, start + window_bp)
            locus_id = f"chr{chrom}_{start + 1:09d}_{end:09d}"
            rows.append(dict(chrom=chrom, start=start, end=end, core_start=start, core_end=end,
                             name=locus_id, locus_id=locus_id,
                             chromosome_window=True))
    return rows


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


def query_rows(vcf: Path, region: str, include_aa: bool = False) -> tuple[list[str], list[list[str]]]:
    samples = list(vcf_samples(vcf))
    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT" + ("\t%INFO/AA" if include_aa else "") + "[\t%GT]\n"
    text = run(["bcftools", "query", "-r", region, "-f", fmt, str(vcf)])
    return samples, [line.split("\t") for line in text.splitlines() if line]


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


def ancestral_allele(value: str, ref: str, alt: str) -> str:
    """Normalize 1KG INFO/AA values such as A, A|, or A|||."""
    token = re.split(r"[|,]", str(value).strip().upper(), maxsplit=1)[0]
    return token if len(token) == 1 and token in BASES and token in {ref, alt} else "N"


def alt_dosage(gt: str, ploidy: int) -> int | None:
    """Return ALT-allele dosage without requiring phased genotypes."""
    if gt in {".", "", "./.", ".|."}:
        return None
    parts = re.split(r"[/|]", gt)
    if len(parts) != ploidy or any(not x.isdigit() or int(x) > 1 for x in parts):
        return None
    return sum(int(x) for x in parts)


def modern_data(
    vcf: Path,
    contig: str,
    locus: dict,
    min_mac: int,
    sexes: dict[str, str],
    x_male_only: bool,
    x_par_diploid: bool,
):
    samples, rows = query_rows(vcf, f"{contig}:{locus['start'] + 1}-{locus['end']}", include_aa=True)
    if not samples:
        fail(f"modern VCF has no samples: {vcf}")
    missing = [sample for sample in samples if sample not in sexes] if sexes else []
    if sexes and missing:
        fail(f"samples absent from samples.txt (first examples): {','.join(missing[:10])}")
    sites = []
    anchors = []
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
        if vid == locus["name"]:
            anchors.append((pos_i, ref, alt))
        if len(ref) != 1 or len(alt) != 1 or ref.upper() not in BASES or alt.upper() not in BASES:
            continue
        haps = []
        for sample, gt in zip(samples, gts):
            if locus["chrom"] == "X":
                if x_male_only:
                    pair = [called_haploid_base(ref.upper(), alt.upper(), gt), "N"]
                elif x_par_diploid:
                    pair = list(called_base(ref.upper(), alt.upper(), gt, phased_required=True))
                else:
                    fail("chrX input requires an explicit male non-PAR or diploid PAR mode")
            else:
                pair = list(called_base(ref.upper(), alt.upper(), gt, phased_required=True))
            haps.extend(pair)
        counts = Counter(x for x in haps if x in BASES)
        if len(counts) < 2 or min(counts.values()) < min_mac:
            continue
        ploidy = 1 if locus["chrom"] == "X" and x_male_only else 2
        dosages = [alt_dosage(gt, ploidy) for gt in gts]
        sites.append(dict(pos=pos_i, vid=vid, ref=ref.upper(), alt=alt.upper(),
                          ancestral=ancestral_allele(aa,ref.upper(),alt.upper()),haps=haps,dosages=dosages))
    if not sites:
        raise SkipLocus(f"no polymorphic biallelic SNPs with MAC >= {min_mac}")
    target_pos = anchors[0][0] if anchors else (locus["core_start"] + locus["core_end"]) // 2
    # LD selection needs an anchor that survived the biallelic/MAC filters.  An
    # exact named anchor is preferred; otherwise use the closest eligible SNP
    # and report its actual position downstream.
    anchor_pos = min(sites, key=lambda x: (abs(x["pos"] - target_pos), x["pos"]))["pos"]
    return samples, sites, anchor_pos


def dosage_r2(anchor: dict, site: dict, n_samples: int, haploid: bool) -> tuple[float | None, int]:
    """Return unphased genotype-dosage r^2, matching PLINK's LD block intent."""
    pairs = [
        (x, y)
        for x, y in zip(anchor["dosages"][:n_samples], site["dosages"][:n_samples])
        if x is not None and y is not None
    ]
    n_complete = len(pairs)
    if n_complete < 3:
        return None, n_complete
    mean_a = sum(x for x, _ in pairs) / n_complete
    mean_b = sum(y for _, y in pairs) / n_complete
    ss_a = sum((x - mean_a) ** 2 for x, _ in pairs)
    ss_b = sum((y - mean_b) ** 2 for _, y in pairs)
    if ss_a <= 0 or ss_b <= 0:
        return None, n_complete
    cross = sum((x - mean_a) * (y - mean_b) for x, y in pairs)
    return max(0.0, min(1.0, (cross * cross) / (ss_a * ss_b))), n_complete


def select_ld_region(
    sites: list[dict],
    anchor_pos: int,
    n_samples: int,
    haploid: bool,
    r2_threshold: float,
    min_ld_sites: int,
) -> tuple[list[dict], dict, list[tuple]]:
    """Select the contiguous span bounded by SNPs in high LD with the anchor."""
    anchor = min(sites, key=lambda x: (abs(x["pos"] - anchor_pos), x["pos"]))
    anchor_pos = anchor["pos"]
    ld_rows = []
    high_ld = []
    for site in sites:
        r2, n_complete = dosage_r2(anchor, site, n_samples, haploid)
        if site["pos"] == anchor_pos:
            r2 = 1.0
        passed = r2 is not None and r2 >= r2_threshold
        if passed:
            high_ld.append(site)
        ld_rows.append((
            site["pos"], site["vid"], site["ref"], site["alt"],
            int(site["pos"] == anchor_pos), n_complete,
            "" if r2 is None else f"{r2:.12g}", int(passed),
        ))
    if len(high_ld) < min_ld_sites:
        raise SkipLocus(
            f"anchor chr position {anchor_pos} has {len(high_ld)} SNP(s) with "
            f"r2 >= {r2_threshold:g}; require at least {min_ld_sites}",
            code="low_ld_information",
            details={
                "anchor_pos": anchor_pos,
                "selection_method": "anchor_unphased_ld_block",
                "ld_r2_threshold": r2_threshold,
                "n_search_sites": len(sites),
                "n_ld_sites": len(high_ld),
                "ld_rows": ld_rows,
            },
        )
    selected_start = min(x["pos"] for x in high_ld) - 1
    selected_end = max(x["pos"] for x in high_ld)
    selected = [x for x in sites if selected_start < x["pos"] <= selected_end]
    meta = {
        "anchor_pos": anchor_pos,
        "selected_start": selected_start,
        "selected_end": selected_end,
        "selected_bp": selected_end - selected_start,
        "selection_method": "anchor_unphased_ld_block",
        "ld_r2_threshold": r2_threshold,
        "n_search_sites": len(sites),
        "n_ld_sites": len(high_ld),
    }
    return selected, meta, ld_rows


def select_full_window(
    sites: list[dict], anchor_pos: int, window_start: int, window_end: int
) -> tuple[list[dict], dict, list[tuple]]:
    """Keep chromosome-scan windows intact; they do not have a flanked BED anchor."""
    ld_rows = [
        (x["pos"], x["vid"], x["ref"], x["alt"], int(x["pos"] == anchor_pos), "", "", 1)
        for x in sites
    ]
    return sites, {
        "anchor_pos": anchor_pos,
        "selected_start": window_start,
        "selected_end": window_end,
        "selected_bp": window_end - window_start,
        "selection_method": "chromosome_window",
        "ld_r2_threshold": None,
        "n_search_sites": len(sites),
        "n_ld_sites": len(sites),
    }, ld_rows


def select_core_region(sites: list[dict], anchor_pos: int, core_start: int, core_end: int) -> tuple[list[dict], dict, list[tuple]]:
    """Use an explicit BED interval as the analysis region without LD shrinkage."""
    selected = [x for x in sites if core_start < x["pos"] <= core_end]
    ld_rows = [
        (x["pos"], x["vid"], x["ref"], x["alt"], int(x["pos"] == anchor_pos), "", "", int(x in selected))
        for x in sites
    ]
    if len(selected) < 2:
        raise SkipLocus(
            f"explicit core interval has {len(selected)} eligible SNP(s); require at least 2",
            code="low_core_information",
            details={
                "anchor_pos": anchor_pos,
                "selected_start": core_start,
                "selected_end": core_end,
                "selected_bp": core_end-core_start,
                "selection_method": "explicit_bed_core",
                "ld_r2_threshold": None,
                "n_search_sites": len(sites),
                "n_ld_sites": len(selected),
                "ld_rows": ld_rows,
            },
        )
    return selected, {
        "anchor_pos": anchor_pos,
        "selected_start": core_start,
        "selected_end": core_end,
        "selected_bp": core_end-core_start,
        "selection_method": "explicit_bed_core",
        "ld_r2_threshold": None,
        "n_search_sites": len(sites),
        "n_ld_sites": len(selected),
    }, ld_rows


def archaic_calls(vcf: Path, contig: str, locus: dict, modern_sites: list[dict]) -> dict[int, str]:
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
        observed = [x for x in observed if x in BASES and x in {wanted[pos_i]["ref"], wanted[pos_i]["alt"]}]
        calls[pos_i] = observed[0] if observed and len(set(observed)) == 1 else "N"
    return calls


def write_tsv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def match_pairs(helper: Path, pair_file: Path, output: Path) -> list[dict]:
    command = 'source "$1"; match_HAP --header "$2"'
    with output.open("w") as handle:
        p = subprocess.run(["bash", "-c", command, "phyml-match", str(helper), str(pair_file)], text=True, stdout=handle)
    if p.returncode:
        fail(f"shared match_HAP failed for {pair_file}")
    with output.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def lineage(ref: str) -> str:
    value = ref.lower()
    if re.search(r"altai|vindija|chagyr|neand", value):
        return "Neanderthal"
    if "denis" in value:
        return "Denisovan"
    return "Archaic"


def analyse_locus(locus: dict, args, refs: list[str]):
    locus_dir = args.out / "loci"
    if not args.flat_locus_dir:
        locus_dir /= locus["locus_id"]
    locus_dir.mkdir(parents=True, exist_ok=True)
    modern_vcf = vcf_path(args.vcf_dir, locus["chrom"])
    modern_contig = vcf_contig(modern_vcf, locus["chrom"])
    samples, search_sites, anchor_pos = modern_data(
        modern_vcf,
        modern_contig,
        locus,
        args.min_mac,
        args.sexes,
        args.x_male_only,
        args.x_par_diploid,
    )
    if locus.get("chromosome_window"):
        sites, selection, ld_rows = select_full_window(
            search_sites, anchor_pos, locus["start"], locus["end"]
        )
    elif args.region_mode == "core":
        try:
            sites, selection, ld_rows = select_core_region(
                search_sites, anchor_pos, locus["core_start"], locus["core_end"]
            )
        except SkipLocus as exc:
            write_tsv(locus_dir / "search_sites.tsv", ["chr", "pos", "id", "ref", "alt"],
                      ((locus["chrom"], x["pos"], x["vid"], x["ref"], x["alt"]) for x in search_sites))
            if exc.details.get("ld_rows") is not None:
                write_tsv(locus_dir / "ld.tsv", ["chr", "pos", "id", "ref", "alt", "is_anchor", "n_complete", "r2", "pass_r2"],
                          ((locus["chrom"], *row) for row in exc.details["ld_rows"]))
            raise
    else:
        try:
            sites, selection, ld_rows = select_ld_region(
                search_sites,
                anchor_pos,
                len(samples),
                locus["chrom"] == "X" and args.x_male_only,
                args.ld_r2,
                args.min_ld_sites,
            )
        except SkipLocus as exc:
            # Preserve the diagnostic evidence even when no defensible tree
            # region can be defined.  A negative-control locus must remain a
            # visible result, not disappear behind a nonzero process exit.
            write_tsv(locus_dir / "search_sites.tsv", ["chr", "pos", "id", "ref", "alt"],
                      ((locus["chrom"], x["pos"], x["vid"], x["ref"], x["alt"]) for x in search_sites))
            if exc.details.get("ld_rows") is not None:
                write_tsv(locus_dir / "ld.tsv", ["chr", "pos", "id", "ref", "alt", "is_anchor", "n_complete", "r2", "pass_r2"],
                          ((locus["chrom"], *row) for row in exc.details["ld_rows"]))
            raise
    anchor_pos = selection["anchor_pos"]
    write_tsv(locus_dir / "search_sites.tsv", ["chr", "pos", "id", "ref", "alt"],
              ((locus["chrom"], x["pos"], x["vid"], x["ref"], x["alt"]) for x in search_sites))
    write_tsv(locus_dir / "ld.tsv", ["chr", "pos", "id", "ref", "alt", "is_anchor", "n_complete", "r2", "pass_r2"],
              ((locus["chrom"], *row) for row in ld_rows))
    selected_locus = {**locus, "start": selection["selected_start"], "end": selection["selected_end"]}
    archaic_sequences = {}
    for ref in refs:
        avcf = archaic_vcf(args.archaic_root, ref, locus["chrom"])
        ac = archaic_calls(avcf, vcf_contig(avcf, locus["chrom"]), selected_locus, sites)
        archaic_sequences[ref] = "".join(ac.get(x["pos"], "N") for x in sites)

    keep = []
    for i, site in enumerate(sites):
        n_called = sum(seq[i] in BASES for seq in archaic_sequences.values())
        if n_called and (not args.require_all_archaic or n_called == len(refs)):
            keep.append(i)
    if not keep:
        raise SkipLocus(
            "no sites callable in the requested archaic references",
            code="low_archaic_information",
            details={
                "anchor_pos": anchor_pos,
                "selected_start": selection["selected_start"],
                "selected_end": selection["selected_end"],
                "selected_bp": selection["selected_bp"],
                "selection_method": selection["selection_method"],
                "ld_r2_threshold": selection["ld_r2_threshold"],
                "n_search_sites": selection["n_search_sites"],
                "n_ld_sites": selection["n_ld_sites"],
            },
        )
    sites = [sites[i] for i in keep]
    archaic_sequences = {ref: "".join(seq[i] for i in keep) for ref, seq in archaic_sequences.items()}
    ancestral_sequence = "".join(x["ancestral"] for x in sites)
    n_ancestral_sites = sum(x in BASES for x in ancestral_sequence)
    selection["n_ancestral_sites"] = n_ancestral_sites
    write_tsv(
        locus_dir / "selected_region.tsv",
        ["locus_id", "chr", "core_start", "core_end", "analysis_start", "analysis_end", "anchor_pos",
         "selected_start", "selected_end", "selected_bp", "selection_method", "ld_r2_threshold",
         "n_search_sites", "n_ld_sites", "n_tree_sites", "n_ancestral_sites"],
        [(
            locus["locus_id"], locus["chrom"], locus["core_start"], locus["core_end"], locus["start"], locus["end"],
            anchor_pos, selection["selected_start"], selection["selected_end"], selection["selected_bp"],
            selection["selection_method"], "NA" if selection["ld_r2_threshold"] is None else selection["ld_r2_threshold"],
            selection["n_search_sites"], selection["n_ld_sites"], len(sites), n_ancestral_sites,
        )],
    )
    print(
        f"PHYML REGION: locus={locus['locus_id']} search={locus['start']}-{locus['end']} "
        f"selected={selection['selected_start']}-{selection['selected_end']} anchor={anchor_pos} "
        f"method={selection['selection_method']} ld_sites={selection['n_ld_sites']} tree_sites={len(sites)}",
        flush=True,
    )

    copies = []
    for sample_i, sample in enumerate(samples):
        # Pure non-PAR male X is haploid (one phased copy per sample).  PAR
        # input is explicitly routed through chrXY and remains diploid.
        hap_indexes = (0,) if locus["chrom"] == "X" and args.x_male_only else (0, 1)
        for hap_i in hap_indexes:
            seq = "".join(site["haps"][2 * sample_i + hap_i] for site in sites)
            copies.append((sample, hap_i + 1, seq))
    grouped: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for sample, hap, seq in copies:
        grouped[seq].append((sample, hap))
    ordered = sorted(grouped.items(), key=lambda item: (-len(item[1]), item[0]))
    hap_rows = []
    for index, (seq, members) in enumerate(ordered, 1):
        hap_id = f"H{index}"
        hap_rows.append(dict(hap_id=hap_id, n=len(members), seq=seq,
                             copies=";".join(f"{s}:{h}" for s, h in members)))

    write_tsv(locus_dir / "sites.tsv", ["chr", "pos", "id", "ref", "alt"],
              ((locus["chrom"], x["pos"], x["vid"], x["ref"], x["alt"]) for x in sites))
    write_tsv(locus_dir / "haplotypes.tsv", ["hap_id", "n", "seq", "copies"],
              ((x["hap_id"], x["n"], x["seq"], x["copies"]) for x in hap_rows))
    write_tsv(locus_dir / "archaic.tsv", ["archaic", "lineage", "seq"],
              ((ref, lineage(ref), seq) for ref, seq in archaic_sequences.items()))
    write_tsv(locus_dir / "ancestral.tsv", ["reference", "n_callable", "seq"],
              [("Ancestral", n_ancestral_sites, ancestral_sequence)])
    pair_file = locus_dir / "match.input.tsv"
    write_tsv(pair_file, ["query_id", "reference_id", "query_haplotype", "reference_haplotype"],
              ((hap["hap_id"], ref, hap["seq"], seq) for hap in hap_rows for ref, seq in archaic_sequences.items()))
    metrics = match_pairs(args.match_helper, pair_file, locus_dir / "match.tsv")
    by_hap = defaultdict(list)
    for metric in metrics:
        by_hap[metric["query_id"]].append(metric)
    result = []
    for hap in hap_rows:
        values = by_hap[hap["hap_id"]]
        best = max(values, key=lambda x: (float(x["prop_match"]) if x["prop_match"] != "NA" else -1,
                                          int(x["n_match"]), x["reference_id"]))
        result.append({**hap, "best_archaic": best["reference_id"], "best_lineage": lineage(best["reference_id"]),
                       "n_compared": int(best["n_compared"]), "n_match": int(best["n_match"]),
                       "prop_match": None if best["prop_match"] == "NA" else float(best["prop_match"])})

    phy = locus_dir / "haplotypes.phy"
    phy_rows = [(x["hap_id"], x["seq"]) for x in result] + [(ref, seq) for ref, seq in archaic_sequences.items()]
    if n_ancestral_sites >= args.min_ancestral_aa:
        phy_rows.append(("Ancestral", ancestral_sequence))
    else:
        print(f"PHYML WARNING: locus={locus['locus_id']} has only {n_ancestral_sites} usable INFO/AA sites; tree remains unrooted",flush=True)
    labels = []
    with phy.open("w") as handle:
        handle.write(f"{len(phy_rows)} {len(sites)}\n")
        for index, (label, seq) in enumerate(phy_rows, 1):
            short = f"S{index:09d}"
            labels.append((short, label))
            handle.write(f"{short:<10} {seq}\n")
    write_tsv(locus_dir / "haplotypes.phy.meta.tsv", ["phy_label", "label"], labels)
    return result, {**locus, **selection, "n_sites": len(sites), "n_haplotypes": len(result)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loci", type=Path)
    parser.add_argument("--loci-map", type=Path)
    parser.add_argument("--chrs", default="")
    parser.add_argument("--vcf-dir", required=True, type=Path)
    parser.add_argument("--archaic-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--match-helper", required=True, type=Path)
    parser.add_argument("--sample-file", type=Path)
    parser.add_argument("--build", choices=["37", "38"], default="37")
    parser.add_argument("--chr-window-bp", type=int, default=500_000)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--refs", default="Altai Chagyr Vindija Denisova")
    parser.add_argument("--min-mac", type=int, default=2)
    parser.add_argument("--min-compared", type=int, default=2)
    parser.add_argument("--min-prop", type=float, default=0.8)
    parser.add_argument("--ld-r2", type=float, default=0.98)
    parser.add_argument("--min-ld-sites", type=int, default=2)
    parser.add_argument("--min-ancestral-aa", type=int, default=2)
    parser.add_argument(
        "--region-mode", choices=("core", "ld"), default="core",
        help="for --loci, analyse the exact BED core (default) or shrink it around the anchor by LD",
    )
    parser.add_argument("--require-all-archaic", action="store_true")
    parser.add_argument(
        "--x-male-only",
        action="store_true",
        help="treat chrX as the pure non-PAR male-only haploid target",
    )
    parser.add_argument(
        "--x-par-diploid",
        action="store_true",
        help="treat chrX loci as diploid PAR data exported from the chrXY pfile",
    )
    parser.add_argument(
        "--flat-locus-dir",
        action="store_true",
        help="write a single resolved locus directly under OUT/loci",
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        fail("--jobs must be positive")
    if not 0 <= args.ld_r2 <= 1:
        fail("--ld-r2 must be between 0 and 1")
    if args.min_ld_sites < 2:
        fail("--min-ld-sites must be at least 2")
    if args.min_ancestral_aa < 1:
        fail("--min-ancestral-aa must be positive")
    if not args.match_helper.is_file():
        fail(f"shared helper not found: {args.match_helper}")
    args.sexes = read_sexes(args.sample_file) if args.sample_file else {}
    if not args.sexes:
        print("PHYML CHECK WARNING: sample/sex metadata absent; chrX.male identity must be guaranteed by GU preflight", flush=True)
    elif any(sex == "unknown" for sex in args.sexes.values()):
        print("PHYML CHECK WARNING: unknown sample sex values exist outside the validated chrX.male target", flush=True)
    refs = [x for x in re.split(r"[, ]+", args.refs) if x]
    if not refs:
        fail("--refs is empty")
    if args.loci:
        loci = read_loci(args.loci, args.loci_map)
    else:
        chrs = [norm_chr(x) for x in re.split(r"[, ]+", args.chrs) if x]
        if not chrs:
            fail("provide --loci or --chrs")
        if args.chr_window_bp < 10_000:
            fail("--chr-window-bp must be at least 10000")
        loci = chromosome_loci(chrs, args.vcf_dir, args.chr_window_bp, args.build)
    if args.x_male_only and args.x_par_diploid:
        fail("--x-male-only and --x-par-diploid are mutually exclusive")
    loci = annotate_loci(loci, args.x_par_diploid)
    if any(locus["chrom"] == "X" for locus in loci) and not (args.x_male_only or args.x_par_diploid):
        fail("chrX requires an explicit pure non-PAR male or diploid PAR target mode")
    if args.flat_locus_dir and len(loci) != 1:
        print(
            "PHYML WARNING: flat per-locus output requires exactly one resolved locus; "
            "using loci/<locus_id>/ because this request resolved to "
            f"{len(loci)} loci/windows",
            flush=True,
        )
        args.flat_locus_dir = False

    unique_chrs = list(dict.fromkeys(locus["chrom"] for locus in loci))
    for chrom in unique_chrs:
        modern = vcf_path(args.vcf_dir, chrom)
        vcf_contig(modern, chrom)
        for ref in refs:
            avcf = archaic_vcf(args.archaic_root, ref, chrom)
            vcf_contig(avcf, chrom)
    if args.check:
        print(f"CHECK PASSED: loci={len(loci)} chromosomes={','.join(dict.fromkeys(x['chrom'] for x in loci))} refs={','.join(refs)}")
        return

    args.out.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out / "input.loci.bed", ["chr", "start", "end", "name"],
              ((x["chrom"], x["start"], x["end"], x["name"]) for x in loci))
    final = args.out / "final"
    final.mkdir(parents=True, exist_ok=True)
    region_header = ["core_start", "core_end", "analysis_start", "analysis_end", "selected_start", "selected_end",
                     "selected_bp", "selection_method", "ld_r2_threshold", "n_search_sites", "n_ld_sites", "n_ancestral_sites"]
    hap_header = ["method", "genome_build", "locus_id", "chr", *region_header,
                  "name", "anchor_pos", "x_class", "hap_id", "n", "best_archaic", "best_lineage", "n_compared", "n_match",
                  "prop_match", "direct_match_pass", "seq", "copies"]
    locus_header = ["method", "genome_build", "locus_id", "chr", *region_header,
                    "name", "anchor_pos", "x_class", "n_sites", "n_haplotypes", "best_archaic", "best_lineage",
                    "n_compared", "n_match", "prop_match", "direct_match_pass", "status"]
    def process(locus):
        try:
            haplotypes, meta = analyse_locus(locus, args, refs)
        except SkipLocus as exc:
            print(f"PHYML SKIP: {locus['locus_id']} status={exc.code}: {exc}", flush=True)
            return locus, None, None, {
                "status": exc.code,
                "reason": str(exc),
                "details": exc.details,
            }
        return locus, haplotypes, meta, None

    executor = ThreadPoolExecutor(max_workers=args.jobs) if args.jobs > 1 else None
    worker = executor.map if executor is not None else map
    n_summaries = n_haplotypes = n_skipped = 0
    # Stream chromosome-window summaries to disk.  A whole chromosome can
    # contain hundreds of windows and many long haplotype strings; retaining
    # every sequence until the last window needlessly scales peak RAM with the
    # chromosome length.  executor.map preserves locus order while this keeps
    # memory bounded by the active worker count.
    with (final / "haplotypes.tsv").open("w", newline="") as hap_handle, \
         (final / "loci.tsv").open("w", newline="") as locus_handle, \
         (final / "haplotype_samples.tsv").open("w", newline="") as sample_handle, \
         (final / "skipped_loci.tsv").open("w", newline="") as skipped_handle:
        hap_writer = csv.writer(hap_handle, delimiter="\t", lineterminator="\n")
        locus_writer = csv.writer(locus_handle, delimiter="\t", lineterminator="\n")
        sample_writer = csv.writer(sample_handle, delimiter="\t", lineterminator="\n")
        skipped_writer = csv.writer(skipped_handle, delimiter="\t", lineterminator="\n")
        hap_writer.writerow(hap_header)
        locus_writer.writerow(locus_header)
        sample_writer.writerow(["method", "genome_build", "locus_id", "hap_id", "sample", "haplotype"])
        skipped_writer.writerow(["locus_id", "chr", "start", "end", "status", "reason"])
        try:
            processed = worker(process, loci)
            for locus, haplotypes, meta, skip_reason in processed:
                if skip_reason is not None:
                    status = skip_reason["status"]
                    detail = skip_reason["details"]
                    skipped_writer.writerow([
                        locus["locus_id"], locus["chrom"], locus["start"], locus["end"],
                        status, skip_reason["reason"],
                    ])
                    locus_writer.writerow([
                        "phyml", f"GRCh{args.build}", locus["locus_id"], locus["chrom"],
                        locus["core_start"], locus["core_end"], locus["start"], locus["end"],
                        detail.get("selected_start", "NA"), detail.get("selected_end", "NA"),
                        detail.get("selected_bp", "NA"), detail.get("selection_method", "NA"),
                        "NA" if detail.get("ld_r2_threshold") is None else detail["ld_r2_threshold"],
                        detail.get("n_search_sites", 0), detail.get("n_ld_sites", 0), detail.get("n_ancestral_sites", 0),
                        locus["name"], detail.get("anchor_pos", "NA"), locus["x_class"],
                        0, 0, "NA", "NA", 0, 0, "NA", 0, status,
                    ])
                    n_skipped += 1
                    continue
                for hap in haplotypes:
                    hap_pass = int(hap["n_compared"] >= args.min_compared and hap["prop_match"] is not None and hap["prop_match"] >= args.min_prop)
                    hap_writer.writerow([
                        "phyml", f"GRCh{args.build}", meta["locus_id"], meta["chrom"], meta["core_start"], meta["core_end"],
                        meta["start"], meta["end"], meta["selected_start"], meta["selected_end"], meta["selected_bp"],
                        meta["selection_method"], "NA" if meta["ld_r2_threshold"] is None else meta["ld_r2_threshold"],
                        meta["n_search_sites"], meta["n_ld_sites"], meta["n_ancestral_sites"], meta["name"], meta["anchor_pos"], meta["x_class"],
                        hap["hap_id"], hap["n"], hap["best_archaic"], hap["best_lineage"], hap["n_compared"], hap["n_match"],
                        "NA" if hap["prop_match"] is None else f"{hap['prop_match']:.12g}", hap_pass, hap["seq"], hap["copies"]
                    ])
                    n_haplotypes += 1
                    for copy in hap["copies"].split(";"):
                        sample, hap_no = copy.rsplit(":", 1)
                        sample_writer.writerow(["phyml", f"GRCh{args.build}", meta["locus_id"], hap["hap_id"], sample, hap_no])
                eligible = [x for x in haplotypes if x["n_compared"] >= args.min_compared and x["prop_match"] is not None]
                best = max(eligible, key=lambda x: (x["prop_match"], x["n_match"], x["n"]), default=None)
                passed = bool(best and best["prop_match"] >= args.min_prop)
                locus_writer.writerow([
                    "phyml", f"GRCh{args.build}", meta["locus_id"], meta["chrom"], meta["core_start"], meta["core_end"],
                    meta["start"], meta["end"], meta["selected_start"], meta["selected_end"], meta["selected_bp"],
                    meta["selection_method"], "NA" if meta["ld_r2_threshold"] is None else meta["ld_r2_threshold"],
                    meta["n_search_sites"], meta["n_ld_sites"], meta["n_ancestral_sites"], meta["name"], meta["anchor_pos"], meta["x_class"],
                    meta["n_sites"], meta["n_haplotypes"],
                    best["best_archaic"] if best else "NA", best["best_lineage"] if best else "NA",
                    best["n_compared"] if best else 0, best["n_match"] if best else 0,
                    f"{best['prop_match']:.12g}" if best else "NA", int(passed),
                    "direct_match_pass" if passed else f"direct_match_below_prop_{args.min_prop}_or_n_{args.min_compared}"
                ])
                n_summaries += 1
        finally:
            if executor is not None:
                executor.shutdown()
    print(f"PHYML DONE: loci={n_summaries} skipped={n_skipped} haplotypes={n_haplotypes} output={args.out}")


if __name__ == "__main__":
    main()
