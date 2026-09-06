#!/usr/bin/env python3
"""Preserve exact named-anchor genotypes, including indels, for GU PhyML loci.

The PhyML alignment intentionally contains single-nucleotide sites only.  This
companion step records the exact BED column-4 marker separately, so an indel
such as ABO rs8176719 is not silently replaced by the nearest eligible SNP.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

# Reuse the tested VCF/locus discovery code from the main module.
from phyml_core import (  # type: ignore
    archaic_vcf,
    lineage,
    query_rows,
    read_loci,
    vcf_contig,
    vcf_path,
)
from comm import enable_wide_csv_fields


enable_wide_csv_fields()


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"ERROR: {message}")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_copies(value: str) -> list[tuple[str, str]]:
    out = []
    for token in str(value).split(";"):
        token = token.strip()
        if not token:
            continue
        sample, hap = token.rsplit(":", 1) if ":" in token else (token, "")
        out.append((sample, hap))
    return out


def locus_dir(out: Path, locus_id: str, flat: bool) -> Path:
    root = out / "loci"
    return root if flat else root / locus_id


def variant_type(ref: str, alts: list[str]) -> str:
    lengths = {len(ref), *(len(alt) for alt in alts)}
    if all(length == 1 for length in lengths):
        return "SNV"
    if len(lengths) == 1:
        return "MNV"
    return "indel"


def anchor_was_provided(locus: dict) -> bool:
    """Distinguish an explicit BED column-4 marker from GU's auto label.

    gu_normalize_loci() creates locus_<chr>_<start>_<end> when column 4 is
    absent.  That label is a locus identifier, not an anchor variant.
    """
    name = str(locus.get("name", ""))
    chrom = str(locus.get("chrom", ""))
    core_start = int(locus.get("core_start", locus.get("start", 0)))
    core_end = int(locus.get("core_end", locus.get("end", 0)))
    expected = f"locus_{chrom}_{core_start}_{core_end}"
    if name == expected:
        return False
    return not bool(re.fullmatch(r"locus_[^_]+_[0-9]+_[0-9]+(?:_[0-9]+)?", name))


def genotype_alleles(gt: str, allele_count: int, haploid: bool) -> list[int | None]:
    gt = str(gt).strip()
    if gt in {"", ".", "./.", ".|."}:
        return [None] if haploid else [None, None]
    if haploid:
        if "/" in gt or "|" in gt:
            parts = re.split(r"[/|]", gt)
            called = [part for part in parts if part not in {"", "."}]
            if len(called) == 2 and called[0] == called[1]:
                parts = [called[0]]
            else:
                return [None]
        else:
            parts = [gt]
    else:
        # Unphased heterozygotes cannot be assigned to haplotype 1 versus 2.
        if "/" in gt:
            parts = gt.split("/")
            if len(parts) == 2 and parts[0] != parts[1]:
                return [None, None]
        else:
            parts = gt.split("|")
        if len(parts) == 1:
            parts *= 2
        if len(parts) != 2:
            return [None, None]
    out: list[int | None] = []
    for part in parts:
        if not part.isdigit() or int(part) >= allele_count:
            out.append(None)
        else:
            out.append(int(part))
    return out


def find_anchor_row(vcf: Path, contig: str, locus: dict) -> tuple[list[str], list[str] | None]:
    samples, rows = query_rows(vcf, f"{contig}:{locus['start'] + 1}-{locus['end']}")
    exact = [row for row in rows if len(row) >= 5 and row[2] == locus["name"]]
    if len(exact) > 1:
        fail(f"{locus['locus_id']}: multiple VCF records have anchor ID {locus['name']}")
    return samples, exact[0] if exact else None


def archaic_anchor_call(vcf: Path, contig: str, pos: int, modern_ref: str, modern_alts: list[str]) -> tuple[str, str, str]:
    """Map an archaic genotype at the exact marker position to modern allele indexes.

    Ancient callsets may contain more than one normalized record at a position.
    Search every record before declaring a representation mismatch.  A
    heterozygous/unphased archaic genotype is deliberately retained as an
    allele set for audit, but it is not converted to one phased haplotype.
    """
    _, rows = query_rows(vcf, f"{contig}:{pos}-{pos}")
    modern_alleles = [modern_ref] + modern_alts
    saw_record = False
    ambiguous_sets: list[str] = []
    mismatches: list[str] = []
    for row in rows:
        if len(row) < 6 or int(row[1]) != pos:
            continue
        saw_record = True
        _, _, _, ref, alt, *gts = row
        alleles = [ref] + alt.split(",")
        observed: list[str] = []
        for gt in gts:
            indexes = genotype_alleles(gt, len(alleles), haploid=False)
            observed.extend(alleles[index] for index in indexes if index is not None)
        called = sorted(set(observed))
        if not called:
            continue
        if ref == modern_ref and all(allele in modern_alleles for allele in called):
            indexes = sorted({modern_alleles.index(allele) for allele in called})
            allele_text = ",".join(called)
            if len(indexes) == 1:
                return str(indexes[0]), allele_text, "mapped"
            ambiguous_sets.append(allele_text)
            continue
        mismatches.append(",".join(called))
    if ambiguous_sets:
        return "", ";".join(sorted(set(ambiguous_sets))), "ambiguous_or_heterozygous"
    if mismatches:
        return "", ";".join(sorted(set(mismatches))), "representation_mismatch"
    return "", "", "not_called" if saw_record else "position_absent"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loci", required=True, type=Path)
    parser.add_argument("--loci-map", type=Path)
    parser.add_argument("--vcf-dir", required=True, type=Path)
    parser.add_argument("--archaic-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--refs", default="Altai Chagyr Vindija Denisova Denisova25")
    parser.add_argument("--common-min-copies", type=int, default=11,
                        help="common-haplotype threshold used by the regional tree (>10 by default)")
    parser.add_argument("--x-male-only", action="store_true")
    parser.add_argument("--x-par-diploid", action="store_true")
    parser.add_argument("--flat-locus-dir", action="store_true")
    args = parser.parse_args()
    if args.common_min_copies < 1:
        fail("--common-min-copies must be positive")

    loci = read_loci(args.loci, args.loci_map)
    refs = [value for value in re.split(r"[, ]+", args.refs) if value]
    summary_rows: list[dict[str, object]] = []
    hap_rows_all: list[dict[str, object]] = []
    copy_rows_all: list[dict[str, object]] = []
    group_rows_all: list[dict[str, object]] = []

    for locus in loci:
        lid = locus["locus_id"]
        ldir = locus_dir(args.out, lid, args.flat_locus_dir)
        provided = anchor_was_provided(locus)
        modern = vcf_path(args.vcf_dir, locus["chrom"])
        if provided:
            samples, anchor = find_anchor_row(modern, vcf_contig(modern, locus["chrom"]), locus)
        else:
            samples, anchor = [], None
        raw_haps = read_tsv(ldir / "haplotypes.tsv")
        raw_hap_by_id = {row.get("hap_id", ""): row for row in raw_haps}
        copy_to_hap: dict[tuple[str, str], str] = {}
        for row in raw_haps:
            for sample, hap in parse_copies(row.get("copies", "")):
                copy_to_hap[(sample, hap)] = row.get("hap_id", "")

        if anchor is None:
            status = "anchor_not_provided" if not provided else "named_anchor_not_found_in_modern_vcf"
            summary_rows.append({
                "locus_id": lid,
                "anchor_provided": int(provided),
                "requested_anchor_id": locus["name"] if provided else "",
                "anchor_id": locus["name"] if provided else "",
                "chr": locus["chrom"],
                "pos": "",
                "ref": "",
                "alt": "",
                "variant_type": "",
                "exact_anchor_found": 0,
                "modern_called_copies": 0,
                "modern_unassigned_copies": 0,
                "archaic_calls": "",
                "neanderthal_consensus_allele_index": "",
                "denisovan_consensus_allele_index": "",
                "status": status,
                "regional_analysis_effect": "none;anchor_is_secondary",
            })
            write_tsv(
                ldir / "anchor.tsv",
                [
                    "locus_id", "anchor_provided", "requested_anchor_id", "anchor_id",
                    "chr", "pos", "ref", "alt", "variant_type", "exact_anchor_found",
                    "modern_called_copies", "modern_unassigned_copies", "archaic_calls",
                    "neanderthal_consensus_allele_index", "denisovan_consensus_allele_index",
                    "status", "regional_analysis_effect",
                ],
                [summary_rows[-1]],
            )
            continue

        _, pos_s, vid, ref, alt, *gts = anchor
        pos = int(pos_s)
        alts = alt.split(",")
        alleles = [ref] + alts
        haploid = locus["chrom"] == "X" and args.x_male_only
        aggregate: dict[str, dict[str, object]] = defaultdict(lambda: {
            "indexes": Counter(), "alleles": Counter(), "called": 0, "unassigned": 0,
        })
        modern_called = modern_unassigned = 0
        for sample, gt in zip(samples, gts):
            indexes = genotype_alleles(gt, len(alleles), haploid=haploid)
            for h_index, allele_index in enumerate(indexes, 1):
                hap_no = "1" if haploid else str(h_index)
                hap_id = copy_to_hap.get((sample, hap_no), "")
                allele = alleles[allele_index] if allele_index is not None else ""
                copy_rows_all.append({
                    "locus_id": lid,
                    "sample": sample,
                    "haplotype": hap_no,
                    "hap_id": hap_id,
                    "anchor_id": vid,
                    "anchor_pos": pos,
                    "allele_index": "" if allele_index is None else allele_index,
                    "allele": allele,
                    "assigned_to_hap_id": int(bool(hap_id)),
                })
                if allele_index is None:
                    continue
                modern_called += 1
                if not hap_id:
                    modern_unassigned += 1
                    continue
                data = aggregate[hap_id]
                data["indexes"][allele_index] += 1  # type: ignore[index]
                data["alleles"][allele] += 1  # type: ignore[index]
                data["called"] = int(data["called"]) + 1

        archaic_details = []
        archaic_by_index: dict[int, list[tuple[str, str]]] = defaultdict(list)
        lineage_indexes: dict[str, list[int]] = defaultdict(list)
        for ref_name in refs:
            avcf = archaic_vcf(args.archaic_root, ref_name, locus["chrom"])
            index_s, allele, status = archaic_anchor_call(
                avcf, vcf_contig(avcf, locus["chrom"]), pos, ref, alts
            )
            lin = lineage(ref_name)
            if index_s != "":
                index_value = int(index_s)
                lineage_indexes[lin].append(index_value)
                archaic_by_index[index_value].append((ref_name, lin))
            archaic_details.append(f"{ref_name}:{index_s or 'NA'}:{allele or 'NA'}:{status}")

        def consensus_index(lin: str, min_refs: int) -> str:
            counts = Counter(lineage_indexes.get(lin, []))
            if not counts:
                return ""
            index, count = counts.most_common(1)[0]
            if count < min_refs or sum(value == count for value in counts.values()) > 1:
                return ""
            return str(index)

        nea_consensus = consensus_index("Neanderthal", 2)
        den_consensus = consensus_index("Denisovan", 1)
        for hap_id, data in sorted(aggregate.items()):
            indexes: Counter = data["indexes"]  # type: ignore[assignment]
            allele_counts: Counter = data["alleles"]  # type: ignore[assignment]
            called_n = int(data["called"])
            index_set = ",".join(str(index) for index in sorted(indexes))
            allele_set = ",".join(sorted(allele_counts))
            mixed = int(len(indexes) > 1)
            raw_hap = raw_hap_by_id.get(hap_id, {})
            hap_n = int(raw_hap.get("n", called_n) or called_n)
            hap_rows_all.append({
                "locus_id": lid,
                "hap_id": hap_id,
                "hap_n": hap_n,
                "copies": raw_hap.get("copies", ""),
                "anchor_id": vid,
                "anchor_pos": pos,
                "ref": ref,
                "alt": alt,
                "variant_type": variant_type(ref, alts),
                "called_copies": called_n,
                "allele_index_set": index_set,
                "allele_set": allele_set,
                "allele_index_counts": ",".join(f"{key}:{value}" for key, value in sorted(indexes.items())),
                "allele_counts": ",".join(f"{key}:{value}" for key, value in sorted(allele_counts.items())),
                "mixed_anchor_within_snp_haplotype": mixed,
                "common_tree_min_copies": args.common_min_copies,
                "common_tree_eligible": int(hap_n >= args.common_min_copies),
                "raw_best_archaic": raw_hap.get("best_archaic", ""),
                "raw_best_lineage": raw_hap.get("best_lineage", ""),
                "raw_n_compared": raw_hap.get("n_compared", ""),
                "raw_n_match": raw_hap.get("n_match", ""),
                "raw_prop_match": raw_hap.get("prop_match", ""),
                "neanderthal_consensus_allele_index": nea_consensus,
                "denisovan_consensus_allele_index": den_consensus,
                "matches_neanderthal_anchor": int(bool(nea_consensus) and not mixed and index_set == nea_consensus),
                "matches_denisovan_anchor": int(bool(den_consensus) and not mixed and index_set == den_consensus),
            })

        # Secondary anchor-stratified zoom.  Every allele is reported; no risk
        # allele is inferred from the BED file, and no allele group can erase
        # the anchor-independent regional inventory.
        locus_copy_rows = [row for row in copy_rows_all if row.get("locus_id") == lid]
        locus_hap_rows = [row for row in hap_rows_all if row.get("locus_id") == lid]
        for allele_index, allele in enumerate(alleles):
            called_rows = [row for row in locus_copy_rows if str(row.get("allele_index", "")) == str(allele_index)]
            assigned_rows = [row for row in called_rows if int(row.get("assigned_to_hap_id", 0)) == 1]
            hap_with_allele = [
                row for row in locus_hap_rows
                if str(allele_index) in str(row.get("allele_index_set", "")).split(",")
            ]
            unambiguous = [row for row in hap_with_allele if int(row.get("mixed_anchor_within_snp_haplotype", 0)) == 0]
            common = [row for row in unambiguous if int(row.get("common_tree_eligible", 0)) == 1]
            archaic_support = archaic_by_index.get(allele_index, [])
            raw_props: list[float] = []
            for row in unambiguous:
                try:
                    value = float(str(row.get("raw_prop_match", "")).strip())
                except (TypeError, ValueError):
                    continue
                if value == value:
                    raw_props.append(value)
            group_rows_all.append({
                "locus_id": lid,
                "anchor_id": vid,
                "anchor_pos": pos,
                "variant_type": variant_type(ref, alts),
                "allele_index": allele_index,
                "allele": allele,
                "n_modern_called_copies": len(called_rows),
                "n_modern_assigned_copies": len(assigned_rows),
                "n_haplotype_types_with_allele": len(hap_with_allele),
                "haplotype_ids_with_allele": ",".join(sorted(str(row.get("hap_id", "")) for row in hap_with_allele)),
                "n_unambiguous_haplotype_types": len(unambiguous),
                "unambiguous_haplotype_ids": ",".join(sorted(str(row.get("hap_id", "")) for row in unambiguous)),
                "n_common_haplotype_types": len(common),
                "common_haplotype_ids": ",".join(sorted(str(row.get("hap_id", "")) for row in common)),
                "n_common_haplotype_copies": sum(int(row.get("hap_n", 0)) for row in common),
                "common_tree_min_copies": args.common_min_copies,
                "mixed_haplotype_types": len(hap_with_allele) - len(unambiguous),
                "raw_best_archaic_references": ",".join(sorted({
                    str(row.get("raw_best_archaic", "")) for row in unambiguous
                    if str(row.get("raw_best_archaic", ""))
                })),
                "raw_best_lineages": ",".join(sorted({
                    str(row.get("raw_best_lineage", "")) for row in unambiguous
                    if str(row.get("raw_best_lineage", ""))
                })),
                "max_raw_prop_match": f"{max(raw_props):.12g}" if raw_props else "",
                "archaic_references_with_allele": ",".join(sorted(ref_name for ref_name, _ in archaic_support)),
                "archaic_lineages_with_allele": ",".join(sorted({lin for _, lin in archaic_support})),
                "interpretation": "anchor_stratified_zoom_only;not_a_regional_discovery_gate",
            })

        summary_rows.append({
            "locus_id": lid,
            "anchor_provided": 1,
            "requested_anchor_id": locus["name"],
            "anchor_id": vid,
            "chr": locus["chrom"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_type": variant_type(ref, alts),
            "exact_anchor_found": 1,
            "modern_called_copies": modern_called,
            "modern_unassigned_copies": modern_unassigned,
            "archaic_calls": ";".join(archaic_details),
            "neanderthal_consensus_allele_index": nea_consensus,
            "denisovan_consensus_allele_index": den_consensus,
            "status": "complete",
            "regional_analysis_effect": "none;anchor_is_secondary",
        })

        write_tsv(
            ldir / "anchor.tsv",
            [
                "locus_id", "anchor_provided", "requested_anchor_id", "anchor_id",
                "chr", "pos", "ref", "alt", "variant_type",
                "exact_anchor_found", "modern_called_copies", "modern_unassigned_copies",
                "archaic_calls", "neanderthal_consensus_allele_index",
                "denisovan_consensus_allele_index", "status", "regional_analysis_effect",
            ],
            [summary_rows[-1]],
        )

    final = args.out / "final"
    write_tsv(
        final / "anchors.tsv",
        [
            "locus_id", "anchor_provided", "requested_anchor_id", "anchor_id",
            "chr", "pos", "ref", "alt", "variant_type",
            "exact_anchor_found", "modern_called_copies", "modern_unassigned_copies",
            "archaic_calls", "neanderthal_consensus_allele_index",
            "denisovan_consensus_allele_index", "status", "regional_analysis_effect",
        ],
        summary_rows,
    )
    write_tsv(
        final / "anchor_haplotypes.tsv",
        [
            "locus_id", "hap_id", "hap_n", "copies", "anchor_id", "anchor_pos", "ref", "alt", "variant_type",
            "called_copies", "allele_index_set", "allele_set",
            "allele_index_counts", "allele_counts",
            "mixed_anchor_within_snp_haplotype", "common_tree_min_copies", "common_tree_eligible",
            "raw_best_archaic", "raw_best_lineage", "raw_n_compared", "raw_n_match", "raw_prop_match",
            "neanderthal_consensus_allele_index",
            "denisovan_consensus_allele_index", "matches_neanderthal_anchor",
            "matches_denisovan_anchor",
        ],
        hap_rows_all,
    )
    write_tsv(
        final / "anchor_allele_groups.tsv",
        [
            "locus_id", "anchor_id", "anchor_pos", "variant_type", "allele_index", "allele",
            "n_modern_called_copies", "n_modern_assigned_copies", "n_haplotype_types_with_allele",
            "haplotype_ids_with_allele", "n_unambiguous_haplotype_types",
            "unambiguous_haplotype_ids", "n_common_haplotype_types", "common_haplotype_ids",
            "n_common_haplotype_copies", "common_tree_min_copies", "mixed_haplotype_types",
            "raw_best_archaic_references", "raw_best_lineages", "max_raw_prop_match",
            "archaic_references_with_allele", "archaic_lineages_with_allele", "interpretation",
        ],
        group_rows_all,
    )
    write_tsv(
        final / "anchor_copies.tsv",
        [
            "locus_id", "sample", "haplotype", "hap_id", "anchor_id", "anchor_pos",
            "allele_index", "allele", "assigned_to_hap_id",
        ],
        copy_rows_all,
    )
    print(
        f"PHYML ANCHOR: loci={len(summary_rows)} exact={sum(int(row['exact_anchor_found']) for row in summary_rows)} "
        f"haplotypes={len(hap_rows_all)} output={final / 'anchors.tsv'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
