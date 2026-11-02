#!/usr/bin/env python3
"""Validate UKB v3 MFI against actual PVAR alleles and export projection INFO.

MFI contains imputation quality and a release-level MAF, not observed allele
counts. No PGEN genotypes, PSAM individuals, or .acount files are read. Inputs
are merged by position within each explicitly declared GRCh37 autosome, never
by row number or rsID. Both input streams must be position sorted unless
--sort-inputs requests a bounded-memory external sort of metadata only.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
import csv
import gzip
import io
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

from analysis_integrity import InputPath, transactional

VERSION = "ukb-mfi-1.0"
BASES = set("ACGT")
MISSING = {"", ".", "NA", "NAN", "NULL", "NONE", "-"}
MFI_FIELDS = ["chrom", "pos", "mfi_line", "mfi_id", "mfi_rsid", "allele1", "allele2",
              "maf_raw", "minor_allele", "info_raw", "input_layout"]
AUDIT_FIELDS = [*MFI_FIELDS, "mfi_records_at_position", "pvar_records_at_position", "pvar_index",
                "pvar_id", "pvar_ref", "pvar_alt", "id_relation", "match_status", "quality_status", "quality_pass"]
MASK_FIELDS = ["chrom", "pos", "pvar_index", "pvar_id", "pvar_ref", "pvar_alt", "mfi_records_at_position",
               "pvar_records_at_position", "mfi_id", "mfi_rsid", "info_raw", "maf_raw", "match_status", "quality_status", "quality_pass"]
SOURCES = {
    "imputation_build_strand": "https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100319",
    "mfi_resource": "https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=1967",
    "v3_file_verified_prefix": "https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_mfi_chr1_v3.txt",
    "legacy_documentation": "https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/ukb_genetic_data_description.txt",
}
LIMITATIONS = [
    "Official UKB v3 MFI data use column 6=MAF and 7=minor allele; the 2018 prose document lists these in the opposite order. Layout is explicit and checked, never silently inferred.",
    "GRCh37 is an explicit input declaration, checked against recognizable PVAR build headers. Absent headers do not verify assembly; no FASTA or liftover check is performed here.",
    "INFO and MAF annotate this imputed release. They do not establish per-genotype quality, observed WGS AC/AN, native joint multiallelic genotypes, or a third allele.",
    "PVAR REF/ALT determine the exported labels; MFI allele1 is not assumed to remain REF after conversion. Only exact or swapped nucleotide pairs are accepted; no strand-complement guess is made.",
    "Every repeated MFI or PVAR position is excluded conservatively, including split multi-ALT records and identical duplicates. This adapter is for unique biallelic SNV projection markers.",
    "target_info.tsv.gz contains quality-passing, unambiguous matches only. Projection must retain its default rejection of unknown INFO; --allow-unknown-info would weaken this mask.",
    "MFI MAF describes the original release and is not recalculated for a PVAR subset or selected participants. Missing INFO/MAF is never replaced with 1 or inferred from counts.",
]


def chromosomes(value):
    result = []
    for part in value.split(","):
        if re.fullmatch(r"\d+-\d+", part):
            lo, hi = map(int, part.split("-")); current = range(lo, hi+1)
            if lo > hi: raise ValueError("Chromosome ranges must increase")
        elif part.isdigit(): current = [int(part)]
        else: raise ValueError("Use autosomes 1-22, e.g. --chromosomes 1-22 or 1,2")
        result.extend(current)
    if not result or len(set(result)) != len(result) or any(c < 1 or c > 22 for c in result):
        raise ValueError("Select unique autosomes 1-22")
    return sorted(result)


def source_paths(raw):
    selected = chromosomes(raw["chromosomes"])
    result = []
    for kind in ("mfi", "pvar"):
        filename, template = raw.get(kind), raw.get(kind+"_template")
        if bool(filename) == bool(template): raise ValueError(f"Supply exactly one --{kind} or --{kind}-template")
        if filename and len(selected) != 1:
            raise ValueError(f"Single --{kind} file needs exactly one --chromosomes value")
        if template and (template.count("{chrom}") != 1 or re.search(r"[{}]", template.replace("{chrom}", ""))):
            raise ValueError(f"--{kind}-template needs exactly one literal {{chrom}} placeholder")
    for chrom in selected:
        row = {"chrom": chrom}
        for kind in ("mfi", "pvar"):
            row[kind] = str(Path(raw.get(kind) or raw[kind+"_template"].format(chrom=chrom)).expanduser().resolve())
        if row["mfi"] == row["pvar"]: raise ValueError("MFI and PVAR cannot be the same file")
        result.append(row)
    return result


def resolve_inputs(raw):
    return [InputPath(row[kind], "file") for row in source_paths(raw) for kind in ("mfi", "pvar")]


@contextmanager
def text_open(path):
    path = Path(path)
    if str(path).endswith(".zst"):
        try: import zstandard
        except ImportError as exc: raise RuntimeError("Reading .zst metadata requires zstandard") from exc
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as stream, io.TextIOWrapper(stream) as handle:
            yield handle
    elif str(path).endswith((".gz", ".bgz")):
        with gzip.open(path, "rt") as handle: yield handle
    else:
        with path.open("r", buffering=8 << 20) as handle: yield handle


def parsed_number(value):
    if value.upper() in MISSING: return None, "missing"
    try: number = float(value)
    except ValueError: return None, "invalid"
    return (number, "valid") if math.isfinite(number) else (None, "missing")


def mfi_rows(path, chrom, layout):
    with text_open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip(): continue
            values = line.split()
            if len(values) != 8: raise ValueError(f"{path}:{line_no}: require exactly 8 headerless MFI columns")
            vid, rsid, raw_pos, a1, a2, sixth, seventh, inf = values
            try: pos = int(raw_pos)
            except ValueError as exc: raise ValueError(f"{path}:{line_no}: invalid MFI position/header") from exc
            if pos < 1: raise ValueError("MFI positions must be 1-based and positive")
            maf, minor = (sixth, seventh) if layout == "v3" else (seventh, sixth)
            value, number_status = parsed_number(maf)
            if number_status == "invalid":
                raise ValueError(f"{path}:{line_no}: invalid MAF for declared --mfi-layout {layout}; v3 columns 6/7 are MAF/minor allele")
            # Some release IDs encode coordinates; when present, they must agree.
            match = re.match(r"^(?:chr)?([0-9]+|X|Y|MT):([0-9]+)(?:[_:]|$)", vid)
            if match and (not match[1].isdigit() or int(match[1]) != chrom or int(match[2]) != pos):
                raise ValueError(f"{path}:{line_no}: MFI encoded ID conflicts with declared chromosome/position")
            yield dict(chrom=str(chrom), pos=pos, mfi_line=line_no, mfi_id=vid, mfi_rsid=rsid,
                       allele1=a1, allele2=a2, maf_raw=maf, minor_allele=minor, info_raw=inf, input_layout=layout)


def check_build_header(line, declared):
    if not re.match(r"##(?:reference|assembly|genome_build|genomeBuild|contig)\s*=", line, re.I): return []
    found = set()
    for build, pattern in (("GRCh37", r"(?<![A-Za-z0-9])(?:GRCh37|hg19|b37|NCBI37)(?!\d)"),
                           ("GRCh38", r"(?<![A-Za-z0-9])(?:GRCh38|hg38|b38|NCBI38)(?!\d)")):
        if re.search(pattern, line, re.I): found.add(build)
    # Bare numeric declarations also occur in converted metadata. Do not infer
    # an assembly from arbitrary contig lengths or numbers in a file path.
    bare = re.match(r"##(?:reference|assembly|genome_build|genomeBuild)\s*=\s*[\"']?(37|38)[\"']?\s*$", line, re.I)
    if bare: found.add("GRCh"+bare[1])
    for match in re.finditer(r"(?:[,<])\s*assembly\s*=\s*[\"']?(37|38)[\"']?(?=[,>])", line, re.I):
        found.add("GRCh"+match[1])
    if any(build != declared for build in found):
        raise ValueError(f"PVAR build header conflicts with declared {declared}: {line.strip()}")
    return sorted(found)


def pvar_rows(path, chrom, declared, audit):
    columns = None; index = 0
    with text_open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip(): continue
            if line.startswith("##"):
                if columns is not None: raise ValueError("PVAR metadata after its column header")
                if re.match(r"##(?:reference|assembly|genome_build|genomeBuild|contig)\s*=", line, re.I):
                    audit["build_header_lines"].append(line.strip())
                    audit["recognized_builds"].extend(check_build_header(line, declared))
                continue
            fields = line.split()
            if columns is None:
                columns = [x.lstrip("#") for x in fields]
                if len(columns) != len(set(columns)) or not {"CHROM", "POS", "ID", "REF", "ALT"} <= set(columns):
                    raise ValueError(f"{path}:{line_no}: PVAR requires unique CHROM POS ID REF ALT columns")
                continue
            if len(fields) != len(columns): raise ValueError(f"{path}:{line_no}: malformed PVAR row")
            raw = dict(zip(columns, fields)); label = raw["CHROM"].removeprefix("chr")
            if label != str(chrom): raise ValueError(f"{path}:{line_no}: PVAR chromosome {label} conflicts with file chromosome {chrom}")
            try: pos = int(raw["POS"])
            except ValueError as exc: raise ValueError(f"{path}:{line_no}: invalid PVAR position") from exc
            if pos < 1: raise ValueError("PVAR positions must be 1-based and positive")
            yield dict(chrom=str(chrom), pos=pos, pvar_index=index, pvar_id=raw["ID"], pvar_ref=raw["REF"], pvar_alt=raw["ALT"])
            index += 1
    if columns is None: raise ValueError(f"Empty PVAR: {path}")


def ordered(rows, label, sorting, directory, memory):
    """Never hold a chromosome-sized metadata table in memory."""
    if sorting:
        unsorted, sorted_path = Path(directory)/(label+".unsorted"), Path(directory)/(label+".sorted")
        with unsorted.open("w", buffering=8 << 20) as handle:
            for row in rows:
                handle.write(str(row["pos"])+"\t"+json.dumps(row, separators=(",", ":"))+"\n")
        subprocess.run(["sort", "--stable", "-t", "\t", "-k1,1n", "-S", f"{memory}M", "-T", str(directory),
                        "-o", str(sorted_path), str(unsorted)], check=True)
        unsorted.unlink()
        with sorted_path.open(buffering=8 << 20) as handle:
            for line in handle: yield json.loads(line.split("\t", 1)[1])
        sorted_path.unlink()
    else:
        previous = 0
        for row in rows:
            if row["pos"] < previous:
                raise ValueError(f"{label} is not position sorted; rerun with --sort-inputs (metadata-only disk sort)")
            previous = row["pos"]
            yield row


def position_groups(rows):
    for pos, group in itertools.groupby(rows, key=lambda r: r["pos"]):
        # Only one coordinate, never a chromosome. A corruption guard also bounds
        # pathological files containing millions of copies of one position.
        values = list(itertools.islice(group, 10001))
        if len(values) > 10000: raise ValueError(f"More than 10000 records at position {pos}; invalid/unsupported metadata")
        yield pos, values


def merged_groups(left, right):
    sentinel = (math.inf, [])
    a, b = next(left, sentinel), next(right, sentinel)
    while a[0] != math.inf or b[0] != math.inf:
        pos = min(a[0], b[0]); aa = a[1] if a[0] == pos else []; bb = b[1] if b[0] == pos else []
        yield pos, aa, bb
        if a[0] == pos: a = next(left, sentinel)
        if b[0] == pos: b = next(right, sentinel)


def compare(mfi, pvar, min_info, min_maf):
    """One conservative decision shared by the MFI audit and PVAR mask."""
    result = dict(mfi_records_at_position=len(mfi), pvar_records_at_position=len(pvar),
                  match_status="unmatched", quality_status="not_evaluated", quality_pass=0)
    if len(mfi) > 1: result["match_status"] = "duplicate_mfi_position"
    elif len(pvar) > 1: result["match_status"] = "duplicate_pvar_position"
    elif not pvar: result["match_status"] = "outside_pvar_subset"
    elif not mfi: result["match_status"] = "mfi_missing"
    else:
        m, p = mfi[0], pvar[0]
        if p["pvar_ref"] not in BASES or p["pvar_alt"] not in BASES or p["pvar_ref"] == p["pvar_alt"]:
            result["match_status"] = "pvar_not_unique_biallelic_snv"
        elif m["allele1"] not in BASES or m["allele2"] not in BASES or m["allele1"] == m["allele2"]:
            result["match_status"] = "mfi_not_unique_biallelic_snv"
        elif {p["pvar_ref"], p["pvar_alt"]} != {m["allele1"], m["allele2"]}:
            result["match_status"] = "allele_pair_mismatch_no_strand_guess"
        else:
            result["match_status"] = "matched_exact" if p["pvar_ref"] == m["allele1"] else "matched_ref_alt_swapped"
            result["id_relation"] = "alternate_id" if p["pvar_id"] == m["mfi_id"] else "rsid" if p["pvar_id"] == m["mfi_rsid"] else "position_and_alleles_only"
            problems = []
            inf, istatus = parsed_number(m["info_raw"]); maf, mstatus = parsed_number(m["maf_raw"])
            if istatus != "valid": problems.append("info_"+istatus)
            elif not 0 <= inf <= 1: problems.append("info_out_of_range")
            elif inf < min_info: problems.append("info_below_threshold")
            if mstatus != "valid": problems.append("maf_"+mstatus)
            elif not 0 <= maf <= .5: problems.append("maf_out_of_range")
            elif maf < min_maf: problems.append("maf_below_threshold")
            if m["minor_allele"].upper() in MISSING: problems.append("minor_allele_missing")
            elif m["minor_allele"] not in {m["allele1"], m["allele2"]}: problems.append("minor_allele_mismatch")
            if p["pvar_id"].upper() in MISSING: problems.append("pvar_id_missing")
            result["quality_status"] = ";".join(problems) if problems else "pass"
            result["quality_pass"] = int(not problems)
    return result


def tsv_writer(handle, fields):
    writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore", restval="NA")
    writer.writeheader()
    return writer


@contextmanager
def info_output(path, tabix):
    if tabix:
        try: import pysam
        except ImportError as exc: raise RuntimeError("--tabix requires pysam") from exc
        with pysam.BGZFile(str(path), "w") as raw, io.TextIOWrapper(raw) as handle: yield handle
    else:
        with gzip.open(path, "wt") as handle: yield handle


@transactional("manifest.json", resolve_inputs=resolve_inputs, input_prefixes=())
def run(args):
    if not math.isfinite(args.min_info) or not 0 <= args.min_info <= 1 or not math.isfinite(args.min_maf) or not 0 <= args.min_maf <= .5:
        raise ValueError("Require min-info in [0,1] and min-maf in [0,.5]")
    if not 16 <= args.sort_memory_mb <= 4096: raise ValueError("sort-memory-mb must be 16..4096")
    paths, out = source_paths(vars(args)), Path(args.out)
    total = Counter(); by_chrom = []; header_audits = []
    info_path = out/"target_info.tsv.gz"
    with gzip.open(out/"mfi_audit.tsv.gz", "wt") as mh, gzip.open(out/"pvar_quality_mask.tsv.gz", "wt") as ph, info_output(info_path, args.tabix) as ih:
        mw, pw, iw = tsv_writer(mh, AUDIT_FIELDS), tsv_writer(ph, MASK_FIELDS), tsv_writer(ih, ["chrom", "pos", "ref", "alt", "info"])
        for source in paths:
            chrom = source["chrom"]; counts = Counter()
            print(f"[dataset] type=UKB_MFI chr={chrom} mfi={source['mfi']} pvar={source['pvar']} layout={args.mfi_layout} genotype_counts=unavailable", flush=True)
            header = dict(chrom=chrom, pvar=source["pvar"], declared_build=args.genome_build,
                          build_header_lines=[], recognized_builds=[])
            with tempfile.TemporaryDirectory(prefix=f"mfi-chr{chrom}-", dir=out) as temporary:
                left = ordered(mfi_rows(source["mfi"], chrom, args.mfi_layout), "mfi", args.sort_inputs, temporary, args.sort_memory_mb)
                right = ordered(pvar_rows(source["pvar"], chrom, args.genome_build, header), "pvar", args.sort_inputs, temporary, args.sort_memory_mb)
                for pos, ms, ps in merged_groups(iter(position_groups(left)), iter(position_groups(right))):
                    decision = compare(ms, ps, args.min_info, args.min_maf)
                    counts["mfi_rows"] += len(ms); counts["pvar_rows"] += len(ps)
                    counts["positions"] += 1
                    counts["position_"+decision["match_status"]] += 1
                    counts["quality_pass"] += decision["quality_pass"]
                    if decision["quality_status"] != "not_evaluated": counts["quality_status_"+decision["quality_status"]] += 1
                    for m in ms:
                        mw.writerow(dict(m, **(ps[0] if len(ps) == 1 else {}), **decision))
                    for p in ps:
                        pw.writerow(dict(p, **(ms[0] if len(ms) == 1 else {}), **decision))
                    if decision["quality_pass"]:
                        iw.writerow(dict(chrom=chrom, pos=pos, ref=ps[0]["pvar_ref"], alt=ps[0]["pvar_alt"], info=ms[0]["info_raw"]))
            if not counts["mfi_rows"] or not counts["pvar_rows"]:
                raise ValueError(f"Chromosome {chrom}: empty MFI or PVAR; no complete quality contract")
            header["recognized_builds"] = sorted(set(header["recognized_builds"]))
            header["verification"] = "declaration_and_compatible_header" if header["recognized_builds"] else "declaration_only_header_uninformative"
            print(f"[build-check] chr={chrom} declared={args.genome_build} header_builds={','.join(header['recognized_builds']) or 'none'} status={header['verification']}", flush=True)
            print(f"[dataset] chr={chrom} mfi_rows={counts['mfi_rows']} pvar_rows={counts['pvar_rows']} quality_pass={counts['quality_pass']} missing_MFI_positions={counts['position_mfi_missing']} INFO_threshold={args.min_info} MAF_threshold={args.min_maf}", flush=True)
            header_audits.append(header); by_chrom.append(dict(chrom=chrom, **counts)); total.update(counts)
    if args.tabix and total["quality_pass"]:
        import pysam
        pysam.tabix_index(str(info_path), seq_col=0, start_col=1, end_col=1, line_skip=1, zerobased=False, force=True)
    fields = ["chrom", *sorted(set().union(*(r.keys() for r in by_chrom)) - {"chrom"})]
    with (out/"summary.tsv").open("w") as handle:
        writer = tsv_writer(handle, fields)
        for row in by_chrom: writer.writerow({k: row.get(k, 0) for k in fields})
        writer.writerow(dict(chrom="all", **{k: total[k] for k in fields if k != "chrom"}))
    return dict(status="complete" if total["quality_pass"] else "no_usable_exact_high_quality_contrast", schema=VERSION,
                genome_build=args.genome_build, mfi_layout=args.mfi_layout, mfi_layout_columns=(
                    "id rsid pos allele1 allele2 maf minor_allele info" if args.mfi_layout == "v3" else
                    "id rsid pos allele1 allele2 minor_allele maf info"),
                counts=dict(total), chromosomes=by_chrom, pvar_build_audits=header_audits,
                source_files=paths, source_documentation=SOURCES, parameters=vars(args), limitations=LIMITATIONS,
                genotype_access="none; PGEN, PSAM and .acount not read", count_evidence="unknown_from_MFI; no observed AC/AN produced",
                target_info_policy="Only unique exact/swapped SNV matches passing all declared quality checks; absent rows must remain unknown/excluded in projection",
                tabix_index="created" if args.tabix and total["quality_pass"] else "not_created_no_pass_rows" if args.tabix else "not_requested",
                memory="Streaming coordinate merge; at most one coordinate group (cap 10000 rows). Optional external sort uses the declared metadata buffer and temporary disk.")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("mfi", "pvar"):
        group = p.add_mutually_exclusive_group(required=True)
        group.add_argument("--"+name, help="Single chromosome metadata file; use one --chromosomes value")
        group.add_argument("--"+name+"-template", help="Path with literal {chrom}; quote it in the shell")
    p.add_argument("--chromosomes", default="1-22")
    p.add_argument("--genome-build", choices=("GRCh37",), required=True, help="Explicit declaration; conflicting PVAR build headers fail")
    p.add_argument("--mfi-layout", choices=("v3", "legacy_doc"), default="v3", help="v3: column6=MAF/7=minor allele; legacy_doc is the opposite order in 2018 prose")
    p.add_argument("--min-info", type=float, default=.8)
    p.add_argument("--min-maf", type=float, default=0., help="Require valid MAF [0,.5]; optionally apply a higher floor")
    p.add_argument("--sort-inputs", action="store_true", help="External bounded-memory sort of metadata if positions are not ordered")
    p.add_argument("--sort-memory-mb", type=int, default=128)
    p.add_argument("--tabix", action="store_true", help="Write BGZF target_info.tsv.gz and its tabix index (requires pysam)")
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
        print(json.dumps(dict(status=result["status"], counts=result["counts"]), sort_keys=True))
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"ukb-mfi: {exc}")


if __name__ == "__main__":
    raise SystemExit(main())
