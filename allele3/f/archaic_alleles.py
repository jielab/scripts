#!/usr/bin/env python3
"""Link predefined nucleotide alleles to the SAME PERSON's GU tract calls.

Consumes existing candidates/targets.tsv or alleles/site_catalog.tsv.gz, GU
native finalized tracts/rosters, and indexed genotype VCFs. No caller is rerun.
IBDmix is unphased: even a phased target GT does not phase the tract. A matching
archaic base and a tract are supporting evidence, never proof of allele origin.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import gzip
import json
import math
from pathlib import Path
import re
import sqlite3
import subprocess
import tempfile

from analysis_integrity import InputPath, identity, transactional
from allele_measurement import phase_relation
from herv_archaic import inspect_dataset, metadata, sample_ids, DEFAULT_REFERENCES
from scan import normalize_file, track_references

BASES = set("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def read_targets(path, maximum):
    result, seen = [], set()
    with text_open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom, pos, ref = row["chrom"].removeprefix("chr"), int(row["pos"]), row["ref"].upper()
            alts = row.get("alts", row.get("alts_by_af", "")).upper().split(",")
            if chrom not in {*(str(i) for i in range(1, 23)), "X", "Y"} or pos < 1:
                raise ValueError("Targets need canonical 1-based SNP coordinates")
            if ref not in BASES or len(alts) < 2 or len(set(alts)) != len(alts) or any(x not in BASES or x == ref for x in alts):
                raise ValueError(f"Expected at least two distinct ALT bases: {chrom}:{pos}")
            if (chrom, pos) in seen:
                raise ValueError(f"Repeated target coordinate: {chrom}:{pos}")
            seen.add((chrom, pos))
            result.append(dict(row, chrom=chrom, pos=pos, ref=ref, alts=alts,
                               site_id=row.get("site_id") or f"{chrom}:{pos}:{ref}"))
            if maximum and len(result) > maximum:
                raise ValueError(f"More than {maximum} target sites; select candidates or explicitly increase --max-targets")
    if not result:
        raise ValueError("No target sites")
    return result


def map_targets(targets, source_build, gu_build, chain, liftover, out):
    """Lift single bases, retain the original identity and complement on '-' only."""
    if source_build == gu_build:
        return [dict(t, gu_chrom=t["chrom"], gu_pos=t["pos"], gu_ref=t["ref"],
                     gu_alts=t["alts"], mapping="same_build", strand="+") for t in targets], []
    if not chain:
        raise ValueError("Different target/GU builds require --target-to-gu-chain (e.g. hg38ToHg19)")
    out = Path(out)
    bed, mapped, missing = out/"targets.source.bed", out/"targets.mapped.bed", out/"targets.unmapped.bed"
    with bed.open("w") as handle:
        for index, t in enumerate(targets):
            handle.write(f"chr{t['chrom']}\t{t['pos']-1}\t{t['pos']}\tt{index}\t0\t+\n")
    subprocess.run([liftover, "-minMatch=1.0", str(bed), str(chain), str(mapped), str(missing)], check=True)
    hits = defaultdict(list)
    with mapped.open() as handle:
        for line in handle:
            chrom, start, end, name, _, strand = line.split()[:6]
            if not re.fullmatch(r"t\d+", name) or int(name[1:]) >= len(targets):
                raise ValueError("Unexpected liftOver target ID")
            hits[int(name[1:])].append((chrom.removeprefix("chr"), int(start), int(end), strand))
    result, audit = [], []
    target_positions = Counter((x[0], x[1]) for group in hits.values() for x in group)
    for index, t in enumerate(targets):
        matches = hits[index]
        reason = "unmapped" if not matches else "ambiguous_mapping" if len(matches) != 1 else "mapped"
        if reason == "mapped":
            chrom, start, end, strand = matches[0]
            if end-start != 1 or chrom != t["chrom"] or strand not in {"+", "-"} or target_positions[chrom, start] != 1:
                reason = "nonunique_or_nonsyntenic_single_base"
            else:
                flip = lambda base: base.translate(COMPLEMENT) if strand == "-" else base
                result.append(dict(t, gu_chrom=chrom, gu_pos=start+1, gu_ref=flip(t["ref"]),
                                   gu_alts=[flip(x) for x in t["alts"]], mapping="unique_single_base_liftover", strand=strand))
        audit.append(dict(site_id=t["site_id"], status=reason))
    return result, audit


def metadata_vcfs(path, modern_template, archaic_template, chrom, references):
    """Read actual source paths from GU provenance; never guess a relocated file."""
    modern, archaic = [], {}
    for line in Path(path).read_text().splitlines():
        fields = line.replace("\\t", "\t").split("\t")
        if not fields:
            continue
        value = fields[-1].rsplit(":", 2)[0]
        if fields[0] == "modern_vcf" or (fields[0] == "modern_source" and len(fields) >= 4 and fields[2] == "vcf"):
            modern.append(value)
        if fields[0] == "archaic_vcf":
            found = [ref for ref in references if ref in Path(value).parts]
            if len(found) == 1:
                archaic[found[0]] = value
    if modern_template:
        modern_path = modern_template.format(chrom=chrom)
    elif len(set(modern)) == 1:
        modern_path = modern[0]
    else:
        raise ValueError("GU modern source is missing or is a disposable PGEN export; provide --modern-vcf-template")
    if archaic_template:
        archaic = {ref: archaic_template.format(chrom=chrom, reference=ref) for ref in references}
    return modern_path, archaic


def decode_base(records, target_ref, base, sample, min_gq=0., min_dp=0.):
    """Measure an explicitly encoded ALT; never invent a missing third allele.

    An ALT's dosage is readable from its own split row, but this does not recover
    a complete multi-allelic genotype. REF is intentionally not reconstructed.
    """
    empty = dict(status="not_encoded", dosage=None, ploidy=None, target_gt_phased=False)
    if not records:
        return dict(empty, status="no_record")
    if any(r.ref != target_ref for r in records):
        return dict(empty, status="ref_mismatch")
    matches = [r for r in records if base in tuple(r.alleles or ())[1:]]
    if len(matches) != 1:
        return dict(empty, status="ambiguous_alt" if matches else "not_encoded")
    record = matches[0]
    alleles = tuple(record.alleles or ())
    if len(alleles) != len(set(alleles)) or any(x not in BASES for x in alleles):
        return dict(empty, status="invalid_alleles")
    if set(record.filter.keys()) - {"PASS", "."}:
        return dict(empty, status="filtered")
    if sample not in record.samples:
        return dict(empty, status="sample_absent")
    call = record.samples[sample]
    gt = call.get("GT")
    if gt is None or len(gt) not in {1, 2} or any(a is None for a in gt):
        return dict(empty, status="missing_gt")
    if any(not isinstance(a, int) or a < 0 or a >= len(alleles) for a in gt):
        return dict(empty, status="invalid_gt")
    if len(records) > 1:
        if any(len(tuple(r.alleles or ())) != 2 for r in records):
            return dict(empty, status="mixed_split_representation")
        if len({r.alleles[1] for r in records}) != len(records):
            return dict(empty, status="duplicated_alt_records")
        split_gts = [r.samples[sample].get("GT") for r in records if sample in r.samples]
        if len(split_gts) == len(records) and all(x is not None and len(x) == len(gt) and all(a in (0, 1) for a in x) for x in split_gts):
            if sum(sum(x) for x in split_gts) > len(gt):
                return dict(empty, status="split_dosage_conflict")
            comparable, _ = phase_relation(records, sample)
            if comparable and any(sum(x[h] for x in split_gts) > 1 for h in range(len(gt))):
                return dict(empty, status="split_phase_conflict")
    for field, limit in (("GQ", min_gq), ("DP", min_dp)):
        if limit <= 0:
            continue
        try:
            value = float(call.get(field))
        except (TypeError, ValueError):
            return dict(empty, status=f"missing_{field.lower()}")
        if not math.isfinite(value) or value < limit:
            return dict(empty, status=f"low_{field.lower()}")
    return dict(status="called", dosage=sum(alleles[a] == base for a in gt), ploidy=len(gt),
                target_gt_phased=bool(getattr(call, "phased", False)),
                representation="split_alt_dosage_only" if len(records) > 1 else "single_record")


def evidence(call, overlap, archaic_call):
    if call["dosage"] is None:
        return "allele_unmeasured"
    if call["dosage"] == 0:
        return "allele_not_carried_in_measured_gt"
    if not overlap:
        return "allele_observed_without_called_tract"
    if archaic_call["dosage"] is None:
        return "same_person_cooccurrence"
    return "cooccurrence_and_archaic_base_match" if archaic_call["dosage"] > 0 else "cooccurrence_without_archaic_base_match"


def membership_database(normalized, targets, connection):
    """A disk-backed, deduplicated target/person/reference relation; no N^2 scan."""
    connection.execute("CREATE TABLE member(pos INTEGER, reference TEXT, sample TEXT, PRIMARY KEY(pos,reference,sample)) WITHOUT ROWID")
    positions = sorted({t["gu_pos"]-1 for t in targets})
    with Path(normalized).open() as handle:
        for line in handle:
            chrom, ref, sample, start, end, *_ = line.rstrip("\n").split("\t")
            lo, hi = bisect_left(positions, int(start)), bisect_left(positions, int(end))
            connection.executemany("INSERT OR IGNORE INTO member VALUES(?,?,?)",
                                   ((pos+1, ref, sample) for pos in positions[lo:hi]))
    connection.commit()


def fetch_records(vcf, chrom, pos):
    contigs = [x for x in vcf.header.contigs if x.removeprefix("chr") == chrom]
    if len(contigs) != 1:
        raise ValueError(f"Expected one chr{chrom} VCF contig")
    return [r for r in vcf.fetch(contigs[0], pos-1, pos) if r.pos == pos]


def resolve_input_files(raw):
    # Templates are paths with fields; protect their parent directories before
    # transaction staging, including native VCFs outside the GU result root.
    result = []
    for key in ("modern_vcf_template", "archaic_vcf_template"):
        if raw.get(key):
            prefix = str(raw[key]).split("{", 1)[0]
            result.append(InputPath(Path(prefix).parent, "directory"))
    root = Path(raw["gu_dataset"])/raw.get("module", "ibdmix")
    for meta in root.glob("chr*/run.meta.tsv"):
        for line in meta.read_text().splitlines():
            fields = line.replace("\\t", "\t").split("\t")
            if fields[0] in {"modern_vcf", "archaic_vcf"} or (fields[0] == "modern_source" and len(fields) >= 4 and fields[2] == "vcf"):
                result.append(InputPath(fields[-1].rsplit(":", 2)[0], "file"))
    return result


@transactional("manifest.json", input_dirs=("gu_dataset",), input_files=("targets", "target_to_gu_chain"), resolve_inputs=resolve_input_files)
def run(args):
    import pysam
    references = track_references(args.references)
    targets = read_targets(args.targets, args.max_targets)
    root, out = Path(args.gu_dataset).resolve(), Path(args.out)
    source = root/args.module
    present_meta = sorted(source.glob("chr*/run.meta.tsv"))
    builds = {metadata(p).get("genome_build", "").removeprefix("GRCh").removeprefix("b") for p in present_meta}
    if len(builds) != 1 or not builds <= {"37", "38"}:
        raise ValueError("GU source must declare one genome build")
    gu_build = builds.pop()
    if args.modern_build != gu_build:
        raise ValueError("--modern-build must match the GU build; target coordinates are mapped, not genotypes")
    mapped, mapping_audit = map_targets(targets, args.target_build, gu_build, args.target_to_gu_chain, args.liftover, out)
    by_chrom = defaultdict(list)
    for target in mapped:
        by_chrom[target["gu_chrom"]].append(target)
    if not by_chrom:
        raise ValueError("No uniquely mapped targets")
    fields = ["site_id", "chrom", "pos", "ref", "allele", "allele_role", "gu_chrom", "gu_pos", "gu_ref", "gu_allele",
              "reference", "sample_id", "tract_call", "genotype_status", "dosage", "ploidy", "target_gt_phased",
              "archaic_gt_status", "archaic_dosage", "evidence", "archaic_haplotype_assignment"]
    summary_fields = ["site_id", "family", "group", "allele", "allele_role", "reference", "tract_call", "n_tested", "n_called",
                      "n_carriers", "AC", "AN", "allele_frequency", "carrier_fraction", "missing_genotypes", "archaic_gt_status", "archaic_dosage"]
    audits, inputs, nrows = [], [identity(args.targets, full=True)], 0
    with gzip.open(out/"individual_evidence.tsv.gz", "wt", newline="") as handle, (out/"allele_tract_summary.tsv").open("w", newline="") as sh:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore"); writer.writeheader()
        sw = csv.DictWriter(sh, summary_fields, delimiter="\t"); sw.writeheader()
        for chrom, items in sorted(by_chrom.items()):
            info = inspect_dataset(root, args.module, (chrom,), references)
            if args.min_lod < info["native_min_lod"] or args.min_length < info["native_min_length"]:
                raise ValueError("Cannot relax finalized GU thresholds")
            people = sorted(sample_ids(info["sample_panel"]))
            meta = source/f"chr{chrom}/run.meta.tsv"
            modern_path, archaic_paths = metadata_vcfs(meta, args.modern_vcf_template, args.archaic_vcf_template, chrom, references)
            inputs += [identity(p) for p in [*info["inputs"], modern_path]]
            with tempfile.TemporaryDirectory(prefix=f"chr{chrom}-", dir=out) as temporary, ExitStack() as stack:
                temporary = Path(temporary); normalized = temporary/"tracts.tsv"
                with normalized.open("w") as norm:
                    for path in info["files"]:
                        normalize_file(path, "ibdmix", norm, args.min_lod, args.min_length, set(people), gu_build, "bed0", (chrom,), references)
                db = sqlite3.connect(temporary/"membership.sqlite"); stack.callback(db.close)
                membership_database(normalized, items, db)
                modern = stack.enter_context(pysam.VariantFile(modern_path))
                archaic = {}
                for ref, path in archaic_paths.items():
                    if Path(path).is_file():
                        av = stack.enter_context(pysam.VariantFile(path))
                        if len(av.header.samples) != 1:
                            raise ValueError(f"Expected one archaic individual in {path}")
                        archaic[ref] = av; inputs.append(identity(path))
                for target in sorted(items, key=lambda t: t["gu_pos"]):
                    records = fetch_records(modern, chrom, target["gu_pos"])
                    members = {ref: set() for ref in references}
                    for ref, person in db.execute("SELECT reference,sample FROM member WHERE pos=?", (target["gu_pos"],)):
                        members[ref].add(person)
                    ancient_records = {ref: fetch_records(vcf, chrom, target["gu_pos"]) for ref, vcf in archaic.items()}
                    for rank, (base, gu_base) in enumerate(zip(target["alts"], target["gu_alts"]), 1):
                        role = f"catalog_alt_{rank}_not_age_rank"
                        calls = {person: decode_base(records, target["gu_ref"], gu_base, person, args.min_gq, args.min_dp) for person in people}
                        for ref in references:
                            ancient_call = decode_base(ancient_records.get(ref, []), target["gu_ref"], gu_base,
                                                       list(archaic[ref].header.samples)[0] if ref in archaic else "", args.min_archaic_gq, args.min_archaic_dp)
                            for inside in (False, True):
                                selected = [p for p in people if (p in members[ref]) == inside]
                                measured = [calls[p] for p in selected if calls[p]["dosage"] is not None]
                                ac, an = sum(x["dosage"] for x in measured), sum(x["ploidy"] for x in measured)
                                sw.writerow(dict(site_id=target["site_id"], family=target.get("family", "."), group=target.get("group", "."),
                                    allele=base, allele_role=role, reference=ref, tract_call="called" if inside else "no_called_tract",
                                    n_tested=len(selected), n_called=len(measured), n_carriers=sum(x["dosage"] > 0 for x in measured), AC=ac, AN=an,
                                    allele_frequency=ac/an if an else "NA", carrier_fraction=sum(x["dosage"] > 0 for x in measured)/len(measured) if measured else "NA",
                                    missing_genotypes=len(selected)-len(measured), archaic_gt_status=ancient_call["status"], archaic_dosage=ancient_call["dosage"] if ancient_call["dosage"] is not None else "NA"))
                            for person in people:
                                call, inside = calls[person], person in members[ref]
                                if not args.emit_all_individuals and not (inside and call["dosage"] is not None and call["dosage"] > 0):
                                    continue
                                writer.writerow(dict(target, allele=base, allele_role=role, gu_allele=gu_base, reference=ref, sample_id=person,
                                    tract_call="called" if inside else "no_called_tract", genotype_status=call["status"],
                                    dosage=call["dosage"] if call["dosage"] is not None else "NA", ploidy=call["ploidy"] if call["ploidy"] is not None else "NA",
                                    target_gt_phased=int(call["target_gt_phased"]), archaic_gt_status=ancient_call["status"],
                                    archaic_dosage=ancient_call["dosage"] if ancient_call["dosage"] is not None else "NA",
                                    evidence=evidence(call, inside, ancient_call), archaic_haplotype_assignment="unavailable_IBDmix_is_unphased"))
                                nrows += 1
                audits.append(dict(chrom=chrom, tested=len(people), missing_from_modern_vcf=len(set(people)-set(modern.header.samples)),
                                   native_references=info["references"], selected_references=list(references), caller_contract=info["caller_contract"],
                                   archaic_vcfs_available=sorted(archaic), n_targets=len(items)))
    (out/"mapping_audit.json").write_text(json.dumps(mapping_audit, indent=2)+"\n")
    if args.target_to_gu_chain:
        inputs.append(identity(args.target_to_gu_chain))
    for before in inputs:
        after = identity(before["path"], full="sha256" in before)
        if after != before:
            raise ValueError("Inputs changed during individual-allele analysis")
    return dict(status="complete", schema="allele3-archaic-alleles-1", n_rows=nrows, chromosome_audits=audits, inputs=inputs,
                settings=vars(args), target_build=args.target_build, gu_build=gu_build,
                individual_table_scope="all_tested_person_allele_reference_rows" if args.emit_all_individuals else "positive_allele_and_tract_cooccurrences_only; summary retains every tested person",
                interpretation="Same-person genotype/tract cooccurrence; base matching supports but cannot establish introgressed origin.",
                limits=["No called tract is not certified absence of ancestry; per-base tract callability is not assumed.",
                        "IBDmix is unphased; target GT phase does not assign an archaic haplotype.",
                        "No selection, causal effect, mutation-age estimate, P value or cross-reference additive ancestry is inferred.",
                        "Within/outside-called-tract frequencies are descriptive and unadjusted for population structure.",
                        "X/Y retain measured haploid/diploid GT ploidy, and are descriptive; caller calibration is not supplied."])


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gu-dataset", required=True, help="Dataset directory, e.g. /mnt/f/analysis/gu/1kg")
    p.add_argument("--module", default="ibdmix", choices=("ibdmix", "ibdmix-GPU"))
    p.add_argument("--targets", required=True, help="Existing candidates/targets.tsv or alleles/site_catalog.tsv.gz")
    p.add_argument("--target-build", choices=("37", "38"), default="38")
    p.add_argument("--modern-build", choices=("37", "38"), required=True, help="Explicit VCF build; must match GU")
    p.add_argument("--target-to-gu-chain", help="For GRCh38 targets + GRCh37 GU, use hg38ToHg19.over.chain.gz")
    p.add_argument("--liftover", default="liftOver")
    p.add_argument("--modern-vcf-template", help="Optional moved VCF override with {chrom}; otherwise GU provenance supplies the path")
    p.add_argument("--archaic-vcf-template", help="Optional override with {reference} and {chrom}; otherwise GU provenance supplies paths")
    p.add_argument("--references", default=",".join(DEFAULT_REFERENCES))
    p.add_argument("--max-targets", type=int, default=5000, help="Candidate follow-up cap; 0 removes it explicitly")
    p.add_argument("--emit-all-individuals", action="store_true", help="Also emit zero-tract/noncarrier/missing individual rows; summary denominators always include them")
    p.add_argument("--min-lod", type=float, default=4.)
    p.add_argument("--min-length", type=int, default=50000)
    for key in ("min-gq", "min-dp", "min-archaic-gq", "min-archaic-dp"):
        p.add_argument("--"+key, type=float, default=0., help="0 means no extra FORMAT quality filter; thresholds >0 require that field")
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.max_targets < 0 or args.min_length < 1 or any(not math.isfinite(x) or x < 0 for x in (args.min_lod, args.min_gq, args.min_dp, args.min_archaic_gq, args.min_archaic_dp)):
        raise ValueError("Invalid threshold")
    print(json.dumps(run(args), indent=2))


if __name__ == "__main__":
    main()
