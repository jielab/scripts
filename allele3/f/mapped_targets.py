#!/usr/bin/env python3
"""Export GRCh37 target labels from an intact completed time-stage mapping.

No liftover, reference lookup, or genotype decoding is repeated. Every requested
GRCh38 nucleotide, including the original REF label, is translated explicitly.
The GRCh37 REF is the reference already checked by the time-stage mapper, which
can differ from the oriented GRCh38 REF. ALT order has no temporal meaning.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path

from analysis_integrity import identity, transactional, validate_manifest

SCHEMA = "allele3-mapped-targets-1"
BASES = frozenset("ACGT")
AUTOSOMES = frozenset(str(c) for c in range(1, 23))
NA = {"", ".", "NA", "N/A", "None", "null"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
TARGET_FIELDS = ["site_id", "chrom", "pos", "ref", "alts", "requested_bases", "build", "family", "group", "match_id",
                 "chrom38", "pos38", "ref38", "requested_alts38", "strand37to38", "reference_changed",
                 "grch37_ref_among_requested_bases"]
BASE_FIELDS = ["site_id", "chrom38", "pos38", "base38", "role38", "chrom37", "pos37", "base37",
               "role37", "strand37to38", "complemented", "reference_changed"]
EXCLUSION_FIELDS = ["site_id", "chrom38", "pos38", "ref38", "requested_alts38", "reason", "detail"]


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else Path(path).open()


def chrom(value):
    value = str(value).removeprefix("chr")
    if value.isdigit():
        value = str(int(value))
    if value not in AUTOSOMES:
        raise ValueError("The completed time-stage mapping contract is autosomal")
    return value


def position(value):
    value = str(value)
    if not value.isdigit() or int(value) < 1:
        raise ValueError("Coordinates must be positive 1-based integers")
    return int(value)


def allele_set(value, ref, label):
    parts = str(value).upper().split(",")
    if not parts or len(parts) != len(set(parts)) or not set(parts) <= BASES or ref in parts:
        raise ValueError(f"Invalid distinct SNV ALT labels in {label}")
    return frozenset(parts)


def fingerprint(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_candidates(path, maximum):
    candidates, seen_coordinates, raw_hashes = {}, {}, {}
    duplicates = 0
    with text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"chrom", "pos", "ref", "alts"} <= set(reader.fieldnames or []):
            raise ValueError("Candidates require chrom,pos,ref,alts and optional site_id")
        for raw in reader:
            c, p, ref = chrom(raw["chrom"]), position(raw["pos"]), raw["ref"].upper()
            if ref not in BASES:
                raise ValueError("Candidate REF must be one canonical nucleotide")
            alts = allele_set(raw["alts"], ref, "candidate")
            site_id = raw.get("site_id", f"{c}:{p}:{ref}")
            if site_id in NA:
                raise ValueError("A supplied candidate site_id cannot be missing")
            for key in ("build", "genome_build", "assembly", "target_build"):
                if raw.get(key) not in NA | {None} and str(raw[key]).lower() not in {"38", "grch38", "hg38"}:
                    raise ValueError("Candidates must be explicitly GRCh38 when a build column is present")
            normalized = dict(raw, site_id=site_id, chrom=c, pos=p, ref=ref, alts=",".join(sorted(alts)))
            digest = fingerprint(normalized)
            if site_id in candidates:
                if raw_hashes[site_id] != digest:
                    raise ValueError("Conflicting duplicate candidate site_id: " + site_id)
                duplicates += 1
                continue
            if (c, p) in seen_coordinates:
                raise ValueError("Multiple candidate site IDs describe the same original coordinate")
            candidates[site_id], raw_hashes[site_id] = normalized, digest
            seen_coordinates[c, p] = site_id
            if len(candidates) > maximum:
                raise ValueError(f"More than {maximum} candidate sites; use a bounded panel or explicitly raise --max-targets")
    if not candidates:
        raise ValueError("No candidate targets")
    return candidates, duplicates


def verified_time(time_dir):
    result = validate_manifest(time_dir, require_evidence=True)
    if not str(result.get("version", "")).startswith("herv-all-allele-time-"):
        raise ValueError("A completed herv_allele_time producer manifest is required")
    if "mapped_catalog.tsv.gz" not in result["output_inventory"]:
        raise ValueError("mapped_catalog.tsv.gz is absent from the verified time output inventory")
    return result


def read_selected_mapping(path, candidates):
    selected, hashes = {}, {}
    counts = Counter()
    with text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"site_id", "chrom", "pos", "ref", "alts", "mapping_status"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Time mapped_catalog lacks original GRCh38 identity/mapping fields")
        for raw in reader:
            counts["catalog_rows_streamed"] += 1
            sid = raw["site_id"]
            if sid not in candidates:
                continue
            counts["selected_catalog_rows"] += 1
            digest = fingerprint(raw)
            if sid in selected:
                if hashes[sid] != digest:
                    raise ValueError("Conflicting duplicate source site_id: " + sid)
                counts["identical_source_duplicates"] += 1
                continue
            selected[sid], hashes[sid] = raw, digest
    return selected, dict(counts)


def oriented(base, strand):
    return base if strand == "+" else base.translate(COMPLEMENT)


def validated_mapping(row):
    """Match herv_time.map_targets, including genuine build-reference changes."""
    required = {"chrom37", "pos37", "ref37", "alts37", "strand37to38"}
    if not required <= set(row):
        raise ValueError("Mapped source row lacks GRCh37 coordinate/allele/strand fields")
    c38, p38, ref38 = chrom(row["chrom"]), position(row["pos"]), row["ref"].upper()
    c37, p37, ref37, strand = chrom(row["chrom37"]), position(row["pos37"]), row["ref37"].upper(), row["strand37to38"]
    if ref38 not in BASES or ref37 not in BASES or strand not in {"+", "-"}:
        raise ValueError("Mapped source row has invalid nucleotide or strand")
    source_alts = allele_set(row["alts"], ref38, "source GRCh38")
    source_bases37 = {oriented(base, strand) for base in {ref38, *source_alts}}
    labels37 = allele_set(row["alts37"], ref37, "source GRCh37")
    if set(labels37) != source_bases37 - {ref37}:
        raise ValueError("Source alts37 disagrees with oriented original REF+ALT labels and the actual ref37")
    changed = int(oriented(ref37, strand) != ref38)
    checks = {"chrom38":c38, "pos38":p38, "ref38":ref38, "source_ref38":ref38,
              "source_to38":"+", "reference_changed":changed}
    for key, expected in checks.items():
        if key in row and row[key] not in NA:
            observed = chrom(row[key]) if key == "chrom38" else position(row[key]) if key == "pos38" else int(row[key]) if key == "reference_changed" else row[key]
            if observed != expected:
                raise ValueError("Source mapping metadata disagrees with its original identity: " + key)
    if row.get("context38") not in NA | {None}:
        context = row["context38"].upper()
        if len(context) != 3 or context[1] != ref38:
            raise ValueError("Source context38 disagrees with ref38")
    return dict(chrom38=c38, pos38=p38, ref38=ref38, chrom37=c37, pos37=p37, ref37=ref37,
                strand37to38=strand, reference_changed=changed, source_bases38={ref38, *source_alts})


def exclusion(target, reason, detail=""):
    return dict(site_id=target["site_id"], chrom38=target["chrom"], pos38=target["pos"], ref38=target["ref"],
                requested_alts38=target["alts"], reason=reason, detail=detail)


def translate(candidates, selected):
    proposed, excluded = {}, []
    for sid, target in candidates.items():
        row = selected.get(sid)
        if row is None:
            excluded.append(exclusion(target, "not_in_time_catalog"))
            continue
        if (chrom(row["chrom"]), position(row["pos"])) != (target["chrom"], target["pos"]):
            excluded.append(exclusion(target, "source38_coordinate_mismatch"))
            continue
        if row["ref"].upper() != target["ref"]:
            excluded.append(exclusion(target, "source38_reference_mismatch", "Source REF="+row["ref"]))
            continue
        source_alts = allele_set(row["alts"], row["ref"].upper(), "source GRCh38")
        requested = {target["ref"], *target["alts"].split(",")}
        missing = requested - {row["ref"].upper(), *source_alts}
        if missing:
            excluded.append(exclusion(target, "requested_bases_missing_from_source", ",".join(sorted(missing))))
            continue
        if row["mapping_status"] != "mapped":
            excluded.append(exclusion(target, "source_mapping_unavailable", row["mapping_status"]))
            continue
        mapping = validated_mapping(row)
        proposed[sid] = dict(mapping, requested_bases38=requested)
    by_coordinate = defaultdict(list)
    for sid, mapping in proposed.items():
        by_coordinate[mapping["chrom37"], mapping["pos37"]].append(sid)
    for coordinate, ids in by_coordinate.items():
        if len(ids) < 2:
            continue
        detail = f"GRCh37 {coordinate[0]}:{coordinate[1]}; candidates=" + ",".join(sorted(ids))
        for sid in ids:
            excluded.append(exclusion(candidates[sid], "mapped_coordinate_collision", detail))
            del proposed[sid]
    targets37, base_rows = [], []
    for sid, mapping in sorted(proposed.items(), key=lambda x: (int(x[1]["chrom37"]), x[1]["pos37"], x[0])):
        target, strand = candidates[sid], mapping["strand37to38"]
        bases37 = {oriented(base, strand) for base in mapping["requested_bases38"]}
        targets37.append(dict(site_id=sid, chrom=mapping["chrom37"], pos=mapping["pos37"], ref=mapping["ref37"],
            alts=",".join(sorted(bases37-{mapping["ref37"]})), requested_bases=",".join(sorted(bases37)), build="GRCh37",
            family=target.get("family"), group=target.get("group"), match_id=target.get("match_id"),
            chrom38=target["chrom"], pos38=target["pos"], ref38=target["ref"], requested_alts38=target["alts"],
            strand37to38=strand, reference_changed=mapping["reference_changed"],
            grch37_ref_among_requested_bases=int(mapping["ref37"] in bases37)))
        for base in [target["ref"], *sorted(mapping["requested_bases38"]-{target["ref"]})]:
            base37 = oriented(base, strand)
            base_rows.append(dict(site_id=sid, chrom38=target["chrom"], pos38=target["pos"], base38=base,
                role38="REF_label" if base == target["ref"] else "ALT", chrom37=mapping["chrom37"], pos37=mapping["pos37"],
                base37=base37, role37="REF_label" if base37 == mapping["ref37"] else "ALT",
                strand37to38=strand, complemented=int(strand == "-"), reference_changed=mapping["reference_changed"]))
    return targets37, base_rows, sorted(excluded, key=lambda x: (int(x["chrom38"]), x["pos38"], x["site_id"]))


def write_table(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:"NA" if row.get(k) is None else row[k] for k in fields})


@transactional("manifest.json", input_dirs=("time_dir",), input_files=("targets",), input_prefixes=())
def run(args):
    if args.max_targets < 1:
        raise ValueError("--max-targets must be a positive explicit limit")
    time_dir = Path(args.time_dir)
    source_manifest = time_dir/"manifest.json"
    mapped = time_dir/"mapped_catalog.tsv.gz"
    print(f"[dataset] mapped-targets: time-dir={time_dir}; candidates={args.targets}", flush=True)
    receipt = verified_time(time_dir)
    print("[build-check] GRCh38 -> GRCh37: reusing the completed time-stage mapping and checked reference labels; no new liftOver, FASTA or PGEN read", flush=True)
    source_before = identity(source_manifest, full=True)
    mapped_before = identity(mapped)
    targets_before = identity(args.targets, full=True)
    candidates, candidate_duplicates = read_candidates(args.targets, args.max_targets)
    selected, scan_counts = read_selected_mapping(mapped, candidates)
    targets37, base_rows, exclusions = translate(candidates, selected)
    if len(targets37) + len(exclusions) != len(candidates):
        raise ValueError("Internal candidate accounting failure")
    if (source_before != identity(source_manifest, full=True) or mapped_before != identity(mapped)
            or targets_before != identity(args.targets, full=True)):
        raise ValueError("Mapping receipt or candidate inputs changed during export")
    out = Path(args.out)
    write_table(out/"targets_grch37.tsv", TARGET_FIELDS, targets37)
    write_table(out/"base_mapping.tsv", BASE_FIELDS, base_rows)
    write_table(out/"exclusions.tsv", EXCLUSION_FIELDS, exclusions)
    status = "not_available" if not targets37 else "complete_with_exclusions" if exclusions else "complete"
    print(f"[mapped-targets] requested={len(candidates)} exported={len(targets37)} excluded={len(exclusions)}; {status}", flush=True)
    print(f"[build-check] mapped-targets: exported={len(targets37)}; excluded={len(exclusions)}; unreported/missing mappings remain unavailable", flush=True)
    return dict(status=status, schema=SCHEMA, n_candidates=len(candidates), n_exported=len(targets37), n_excluded=len(exclusions),
                n_base_mappings=len(base_rows), identical_candidate_duplicates=candidate_duplicates, stream_counts=scan_counts,
                exclusion_reasons=dict(Counter(r["reason"] for r in exclusions)), source_time_status=receipt["status"],
                source_time_version=receipt["version"], source_time_manifest=source_before, source_mapped_catalog=mapped_before,
                candidate_identity=targets_before, source_build="GRCh38", target_build="GRCh37",
                source_validation="Completed time output_inventory checked before reading; mapped_catalog must be inventoried; small manifest/candidate SHA and mapped-catalog size/mtime rechecked after streaming",
                selected_record_validation="Original coordinates and REF must match; every requested nucleotide must be in the source; selected duplicate IDs must have identical records",
                mapping_semantics="Reuse time mapper's reciprocal single-base mapping and FASTA37 REF; orient every original REF/ALT nucleotide; no new liftOver or reference query",
                requested_bases_contract="The requested_bases column is the exact oriented original REF+requested ALT set; a structural GRCh37 REF outside it must not become an additional measurement",
                limits=["A mapping/reference label does not establish that an allele is observed or encoded in the UKB panel.",
                        "Original GRCh38 site_id is preserved; read chrom/pos/ref fields for GRCh37 identity.",
                        "ALT output is a sorted set, never mutation order, age, priority, or allele-frequency rank.",
                        "The GRCh37 REF can be outside the requested modern bases; it is an assembly label, not a newly observed allele.",
                        "Unmapped and colliding candidates are excluded, never converted into zero counts.",
                        "Only selected candidate records receive semantic checks; unrelated catalog rows are streamed and discarded."])


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--time-dir", required=True, help="Completed time stage containing its verified manifest.json and mapped_catalog.tsv.gz")
    p.add_argument("--targets", required=True, help="Existing GRCh38 candidates/targets.tsv")
    p.add_argument("--max-targets", type=int, default=5000)
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    result = run(parser().parse_args(argv))
    return 0 if result["n_exported"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
