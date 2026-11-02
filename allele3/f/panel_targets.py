#!/usr/bin/env python3
"""Count prespecified nucleotide alleles in an autosomal PGEN candidate panel.

Targets require exact chrom/pos/ref/alts and the same declared assembly as the
panel. This is bounded candidate follow-up, not a variant discovery scan. Every
same-position PVAR companion is retained before applying quality. Counts use
PGEN hardcalls and actual PSAM sample numbers; stored posterior dosages are not
read. A missing/unencoded allele is unavailable, never imputed as zero.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
import os
from pathlib import Path
import re

import numpy as np
from allele_measurement import BASES, PanelRecord, dosage_plan, read_counted_dosage, panel_pvar
from analysis_integrity import InputPath, identity, transactional, validate_manifest
from ukb_mfi import check_build_header, text_open

VERSION = "panel-targets-1.0"
MISSING = {"", ".", "NA", "nan", "None"}
FIELDS = ["site_id", "chrom", "pos", "target_build", "target_ref", "target_alts", "requested_bases", "structural_ref_requested", "allele", "allele_role",
          "population", "group_kind", "evidence_kind", "measurement_status", "counting_scope", "pvar_index",
          "panel_alleles", "allele_idx", "reference_label_relation", "n_panel_samples", "n_selected_samples",
          "n_called", "n_missing", "n_carriers", "AC", "AN", "AF", "call_rate"]
AVAILABILITY = ["site_id", "chrom", "pos", "target_ref", "target_alts", "requested_bases", "structural_ref_requested", "panel_record_count", "panel_representation",
                "invalid_reason", "reference_label_relation", "n_requested_bases", "n_measurable_bases",
                "measurable_bases", "unavailable_bases", "joint_multiallelic_inference"]
LIMITS = [
    "Counts are measured PGEN hardcalls from the declared evidence source. Imputed hardcalls are not posterior dosages or observed WGS chromosomes; this module never calls read_dosages.",
    "Every base has its own callable AN. Missing encoded genotypes and unencoded bases are distinct; sums of split-ALT counts do not establish a native joint multiallelic genotype.",
    "A single biallelic row may be an incomplete split. Its REF is not reconstructed as two minus ALT. REF is counted only when directly encoded as an ALT or in an eligible native multiallelic row.",
    "When targets declare requested_bases, only that exact subset is counted. A structural assembly REF outside this subset is retained for representation checks but never added as a requested or observed allele.",
    "All same-position PVAR companions are checked before quality filtering. A hidden low-quality or non-SNV companion cannot create a false native-REF measurement.",
    "No liftover, strand complement, mutation-age inference, selection test, or new third-allele discovery is performed.",
    "Autosomes only; diploid hardcall AN=2*n_called. X/Y are rejected explicitly because their actual sample ploidy requires a separate contract.",
    "Population labels are supplied by the user. ALL is the selected pooled cohort, not ancestry-adjusted; stratified frequencies do not remove relatedness, ascertainment or batch effects.",
    "Build declarations and recognized PVAR headers are checked, but a FASTA reference audit belongs to the upstream data audit.",
]


def build(value):
    value = value.removeprefix("GRCh")
    if value not in {"37", "38"}: raise ValueError("Require explicit GRCh37 or GRCh38")
    return "GRCh"+value


def chrom(value):
    value = str(value).removeprefix("chr")
    if value not in {str(i) for i in range(1, 23)}:
        raise ValueError("panel targets support autosomes 1-22 only; X/Y require actual ploidy and are explicitly unsupported")
    return value


def table(path, required):
    with text_open(path) as handle:
        names = None
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("##"): continue
            values = line.split()
            if names is None:
                names = [v.lstrip("#") for v in values]
                if len(set(names)) != len(names) or not set(required) <= set(names):
                    raise ValueError(f"{path}: need unique columns {','.join(required)}")
                continue
            if len(values) != len(names): raise ValueError(f"{path}:{line_no}: malformed table row")
            yield dict(zip(names, values))
        if names is None: raise ValueError(f"Empty table: {path}")


def load_targets(path, maximum, declared_build=None):
    result, ids, coordinates = [], set(), set()
    for raw in table(path, ("site_id", "chrom", "pos", "ref", "alts")):
        for column in ("build", "genome_build", "assembly", "target_build"):
            if raw.get(column) not in MISSING | {None}:
                observed = build(raw[column])
                if declared_build is not None and observed != build(declared_build):
                    raise ValueError(f"Target {raw['site_id']} {column}={raw[column]} conflicts with --target-build {declared_build}")
        c, pos, ref, alts = chrom(raw["chrom"]), int(raw["pos"]), raw["ref"].upper(), raw["alts"].upper().split(",")
        sid = raw["site_id"]
        if sid in MISSING or sid in ids or (c, pos) in coordinates:
            raise ValueError("Targets need unique nonmissing site_id and unique coordinates")
        if pos < 1 or ref not in BASES or not 1 <= len(alts) <= 3 or len(set(alts)) != len(alts) or any(a not in BASES or a == ref for a in alts):
            raise ValueError("Targets require 1-based coordinates and distinct A/C/G/T REF/ALTs")
        requested = raw["requested_bases"].upper().split(",") if "requested_bases" in raw else [ref, *alts]
        if not requested or len(requested) != len(set(requested)) or not set(requested) <= {ref, *alts}:
            raise ValueError("requested_bases must be a nonempty distinct subset of the structural REF+ALT labels")
        ids.add(sid); coordinates.add((c, pos))
        result.append(dict(site_id=sid, chrom=c, pos=pos, ref=ref, alts=alts, requested_bases=requested))
        if maximum and len(result) > maximum: raise ValueError(f"Target cap {maximum} exceeded; choose candidates or explicitly increase --max-targets")
    if not result: raise ValueError("No target sites")
    return sorted(result, key=lambda r: (int(r["chrom"]), r["pos"]))


def template(value, chromosome, name):
    if value.count("{chrom}") != 1 or re.search(r"[{}]", value.replace("{chrom}", "")):
        raise ValueError(f"{name} must contain one literal {{chrom}} placeholder")
    return Path(value.format(chrom=chromosome)).expanduser().resolve()


def source_paths(raw):
    targets = load_targets(raw["targets"], raw["max_targets"], raw["target_build"])
    sources = []
    for c in sorted({t["chrom"] for t in targets}, key=int):
        prefix = template(raw["pgen_template"], c, "--pgen-template")
        if prefix.suffix == ".pgen": raise ValueError("--pgen-template is a prefix without .pgen")
        pvar = template(raw["pvar_template"], c, "--pvar-template") if raw.get("pvar_template") else panel_pvar(prefix)
        sources.append(dict(chrom=c, pgen=str(prefix)+".pgen", pvar=str(pvar), psam=str(prefix)+".psam"))
    return targets, sources


def resolve_inputs(raw):
    _, sources = source_paths(raw)
    result = [InputPath(row[k], "file") for row in sources for k in ("pgen", "pvar", "psam")]
    if raw.get("quality_info"):
        parent = Path(raw["quality_info"]).expanduser().resolve().parent
        if (parent/"manifest.json").is_file(): result.append(InputPath(parent, "directory"))
    return result


def quality_contract(args, sources):
    if not args.quality_info:
        if args.evidence_kind != "direct_gt": raise ValueError("imputed_hardcall requires --quality-info; PGEN alone is not imputation QC")
        return dict(status="direct_gt_declared_no_INFO", interpretation="User-declared direct genotype source; PVAR FILTER checked, no per-genotype read-depth certification")
    quality = Path(args.quality_info).expanduser().resolve(); manifest = quality.parent/"manifest.json"
    if manifest.is_file():
        proof = validate_manifest(quality.parent, require_evidence=True)
        if proof.get("schema") != "ukb-mfi-1.0" or quality.name != "target_info.tsv.gz":
            raise ValueError("Existing quality manifest is not the UKB MFI target_info contract")
        if proof.get("genome_build") != build(args.panel_build): raise ValueError("Quality INFO build differs from declared panel build")
        for source in sources:
            path = str(Path(source["pvar"]).resolve()); saved = proof.get("inputs", {}).get(path)
            if not saved: raise ValueError("MFI quality does not bind this exact PVAR: "+path)
            current = identity(path, full="sha256" in saved)
            keys = ("size", "sha256") if "sha256" in saved else ("size", "mtime_ns")
            if any(current.get(k) != saved.get(k) for k in keys): raise ValueError("PVAR changed since MFI conversion: "+path)
        if args.quality_build and build(args.quality_build) != proof["genome_build"]:
            raise ValueError("Explicit --quality-build conflicts with MFI manifest")
        return dict(status="verified_mfi_output_and_original_PVAR", manifest=identity(manifest, full=True), genome_build=proof["genome_build"])
    if not args.quality_build: raise ValueError("Standalone INFO without an MFI manifest requires explicit --quality-build")
    if build(args.quality_build) != build(args.panel_build): raise ValueError("Quality INFO build differs from declared panel build")
    return dict(status="standalone_INFO_build_declared_not_source_verified", genome_build=build(args.quality_build))


def load_quality(path, targets):
    if not path: return {}
    wanted = {(t["chrom"], t["pos"]) for t in targets}; result = {}
    for row in table(path, ("chrom", "pos", "ref", "alt", "info")):
        coordinate = chrom(row["chrom"]), int(row["pos"])
        if coordinate not in wanted: continue
        ref, alt = row["ref"], row["alt"]
        if ref not in BASES or alt not in BASES or ref == alt: raise ValueError("Quality table requires individual SNV contrasts, not a comma-joined ALT")
        key = (*coordinate, ref, alt)
        if key in result: raise ValueError("Duplicate quality contrast: "+str(key))
        value = None if row["info"] in MISSING else float(row["info"])
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1): raise ValueError("INFO must be missing or finite in [0,1]")
        result[key] = value
    return result


def index_pvar(path, chromosome, targets, quality, use_quality, min_info, declared):
    """One metadata pass: exact row indexes plus compact native allele offsets."""
    wanted = {t["pos"] for t in targets}; found = {p: [] for p in wanted}; invalid = {}; nrecords = Counter()
    offsets = array("Q", [0]); names = None; n = 0; headers = []; header_builds = set()
    with text_open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip(): continue
            if line.startswith("##"):
                if names is not None: raise ValueError("PVAR metadata appears after its column header")
                builds = check_build_header(line, declared); header_builds.update(builds)
                if builds: headers.append(line.strip())
                continue
            fields = line.split()
            if names is None:
                names = [v.lstrip("#") for v in fields]
                if len(names) != len(set(names)) or not {"CHROM", "POS", "ID", "REF", "ALT"} <= set(names):
                    raise ValueError("PVAR requires unique CHROM POS ID REF ALT columns")
                continue
            if len(fields) != len(names): raise ValueError(f"{path}:{line_no}: malformed PVAR row")
            row = dict(zip(names, fields)); c, pos = chrom(row["CHROM"]), int(row["POS"])
            if c != chromosome or pos < 1: raise ValueError("PVAR chromosome/position conflicts with its declared template chromosome")
            alleles = (row["REF"], *row["ALT"].split(","))
            if any(not a for a in alleles): raise ValueError("Empty PVAR allele label makes allele offsets unreliable")
            offsets.append(offsets[-1]+len(alleles))
            if pos in wanted:
                nrecords[pos] += 1
                if nrecords[pos] > 64:
                    invalid[pos] = "too_many_same_position_companions"
                elif len(set(alleles)) != len(alleles) or not set(alleles) <= BASES or not 2 <= len(alleles) <= 4:
                    invalid[pos] = "non_snv_or_invalid_companion_at_locus"
                else:
                    hq = frozenset(a for a in alleles[1:] if not use_quality or (
                        quality.get((c, pos, alleles[0], a)) is not None and quality[c, pos, alleles[0], a] >= min_info))
                    found[pos].append(PanelRecord(c, pos, n, alleles, row.get("FILTER", "PASS") in {"PASS", "."}, hq, row["ID"]))
            n += 1
    if names is None or not n: raise ValueError("PVAR has no variants")
    return found, invalid, nrecords, np.frombuffer(offsets, dtype=np.uint64), dict(
        n_pvar_rows=n, allele_offset_bytes=len(offsets)*offsets.itemsize, declared_build=declared,
        build_headers=headers, recognized_builds=sorted(header_builds),
        build_verification="declaration_and_compatible_header" if header_builds else "declaration_only_header_uninformative")


def load_sample_selection(path):
    if not path: return None
    result = {}
    for row in table(path, ("sample_id",)):
        sid = row["sample_id"]
        if sid in MISSING or sid in result: raise ValueError("Samples TSV needs unique sample_id values (PSAM IID)")
        population = row.get("population")
        if population is not None and population in MISSING: population = "UNKNOWN"
        result[sid] = population
    if not result: raise ValueError("Empty sample selection")
    return result


def sample_groups(psam, selection):
    ids = [r["IID"] for r in table(psam, ("IID",))]
    if not ids or len(set(ids)) != len(ids) or any(i in MISSING for i in ids):
        raise ValueError("PSAM needs nonmissing unique IID; ambiguous FID/IID pairs cannot be matched by sample_id")
    lookup = {sid: i for i, sid in enumerate(ids)}
    if selection is not None and not set(selection) <= lookup.keys():
        raise ValueError("Selected sample IDs are absent from PSAM")
    chosen = ids if selection is None else [sid for sid in ids if sid in selection]
    groups = [("ALL", "selected_pooled_unadjusted", np.asarray([lookup[s] for s in chosen], dtype=np.int64))]
    if selection is not None and any(v is not None for v in selection.values()):
        populations = defaultdict(list)
        for sid in chosen: populations[selection[sid] or "UNKNOWN"].append(lookup[sid])
        groups.extend((pop, "user_population_label", np.asarray(indices, dtype=np.int64)) for pop, indices in sorted(populations.items()))
    return len(ids), groups


def plan_target(target, records, invalid_reason, nrecords):
    native = len(records) == 1 and len(records[0].alleles) > 2 and nrecords == 1
    representation = "native_multiallelic_row" if native else "biallelic_or_split_rows" if records and all(len(r.alleles) == 2 for r in records) else "mixed_or_unavailable"
    panel_refs = {r.alleles[0] for r in records}; all_bases = set().union(*(set(r.alleles) for r in records)) if records else set()
    relation = "same_REF_label" if panel_refs == {target["ref"]} else "target_REF_explicitly_in_panel_alleles" if target["ref"] in all_bases else "target_REF_not_encoded"
    reason = invalid_reason or ("target_reference_not_represented" if records and target["ref"] not in all_bases else None)
    requested = set(target.get("requested_bases", [target["ref"], *target["alts"]]))
    entries = []
    for index, base in enumerate([target["ref"], *target["alts"]]):
        if base not in requested: continue
        plan = dict(status=reason) if reason else dosage_plan(records, base)
        entries.append(dict(site_id=target["site_id"], chrom=target["chrom"], pos=target["pos"],
            target_ref=target["ref"], target_alts=",".join(target["alts"]), allele=base,
            requested_bases=",".join(sorted(requested)), structural_ref_requested=int(target["ref"] in requested),
            allele_role="catalog_REF" if index == 0 else f"catalog_ALT_{index}_not_age_rank",
            measurement_plan=plan, reference_label_relation=relation))
    available = [e["allele"] for e in entries if e["measurement_plan"]["status"] == "assayed_pass"]
    audit = dict(site_id=target["site_id"], chrom=target["chrom"], pos=target["pos"], target_ref=target["ref"], target_alts=",".join(target["alts"]),
                 requested_bases=",".join(sorted(requested)), structural_ref_requested=int(target["ref"] in requested),
                 panel_record_count=nrecords, panel_representation=representation, invalid_reason=reason or "none", reference_label_relation=relation,
                 n_requested_bases=len(entries), n_measurable_bases=len(available), measurable_bases=",".join(available) or "NA",
                 unavailable_bases=",".join(e["allele"] for e in entries if e["allele"] not in available) or "NA",
                 joint_multiallelic_inference="native_encoding_only_not_new_variant_discovery" if native and not reason else "not_established_by_separate_ALT_counts")
    return entries, audit


def count_rows(entry, counts, groups, total_samples, args):
    plan = entry["measurement_plan"]; available = plan["status"] == "assayed_pass"
    base = dict(entry, target_build=build(args.target_build), evidence_kind=args.evidence_kind,
                measurement_status=plan["status"], counting_scope=plan.get("measurement_scope", "unavailable"),
                pvar_index=plan.get("pvar_index", "NA"), allele_idx=plan.get("allele_idx", "NA"),
                panel_alleles=",".join(plan.get("panel_alleles", [])) or "NA", n_panel_samples=total_samples)
    for population, kind, indices in groups:
        row = dict(base, population=population, group_kind=kind, n_selected_samples=len(indices))
        if available:
            values = counts[indices]; called = values >= 0; n_called = int(called.sum())
            ac = int(values[called].sum()); an = 2*n_called
            row.update(n_called=n_called, n_missing=len(indices)-n_called, n_carriers=int(np.sum(values[called] > 0)),
                       AC=ac, AN=an, AF=ac/an if an else "NA", call_rate=n_called/len(indices))
        else:
            row.update({k: "NA" for k in ("n_called", "n_missing", "n_carriers", "AC", "AN", "AF", "call_rate")})
        yield row


def writer(handle, fields):
    result = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore", restval="NA")
    result.writeheader(); return result


@transactional("manifest.json", input_files=("targets", "samples", "quality_info"), resolve_inputs=resolve_inputs, input_prefixes=())
def run(args):
    if build(args.target_build) != build(args.panel_build): raise ValueError("Targets and panel builds must match; use an independently audited mapping before this module")
    if args.max_targets < 0 or not 1 <= args.chunk_size <= 256 or not math.isfinite(args.min_info) or not 0 <= args.min_info <= 1:
        raise ValueError("Invalid target cap/chunk size/minimum INFO")
    targets, sources = source_paths(vars(args)); quality_proof = quality_contract(args, sources)
    quality = load_quality(args.quality_info, targets); selection = load_sample_selection(args.samples)
    by_chrom = defaultdict(list)
    for t in targets: by_chrom[t["chrom"]].append(t)
    out = Path(args.out); total = Counter(); chromosome_audits = []; saved_plans = []
    with gzip.open(out/"allele_counts.tsv.gz", "wt") as ch, (out/"target_availability.tsv").open("w") as ah, (out/"sample_groups.tsv").open("w") as sh:
        cw, aw, sw = writer(ch, FIELDS), writer(ah, AVAILABILITY), writer(sh, ["chrom", "population", "group_kind", "n_panel_samples", "n_selected_samples"])
        for source in sources:
            c = source["chrom"]; n_samples, groups = sample_groups(source["psam"], selection)
            print(f"[dataset] chr={c} pgen={source['pgen']} pvar={source['pvar']} psam={source['psam']} evidence={args.evidence_kind} actual_PSAM_N={n_samples} selected_N={len(groups[0][2])}", flush=True)
            for pop, kind, indices in groups:
                sw.writerow(dict(chrom=c, population=pop, group_kind=kind, n_panel_samples=n_samples, n_selected_samples=len(indices)))
            found, invalid, nrecords, offsets, audit = index_pvar(source["pvar"], c, by_chrom[c], quality,
                bool(args.quality_info), args.min_info, build(args.panel_build))
            print(f"[build-check] chr={c} targets={build(args.target_build)} panel={build(args.panel_build)} header_builds={','.join(audit['recognized_builds']) or 'none'} status={audit['build_verification']} quality={quality_proof['status']}", flush=True)
            selected = []; reader_audit = Counter()
            for target in by_chrom[c]:
                entries, availability = plan_target(target, found[target["pos"]], invalid.get(target["pos"]), nrecords[target["pos"]])
                aw.writerow(availability); total["target_sites"] += 1
                for entry in entries:
                    saved_plans.append(entry); total["requested_nucleotides"] += 1
                    total["status_"+entry["measurement_plan"]["status"]] += 1
                    if entry["measurement_plan"]["status"] == "assayed_pass": selected.append(entry)
                    else:
                        for row in count_rows(entry, None, groups, n_samples, args): cw.writerow(row)
            selected.sort(key=lambda e: (e["measurement_plan"]["pvar_index"], e["measurement_plan"]["allele_idx"]))
            if selected:
                import pgenlib
                with pgenlib.PgenReader(os.fsencode(source["pgen"]), raw_sample_ct=n_samples, allele_idx_offsets=offsets) as reader:
                    if reader.get_raw_sample_ct() != n_samples or reader.get_variant_ct() != len(offsets)-1:
                        raise ValueError("PGEN dimensions differ from the actual PSAM/PVAR")
                    for start in range(0, len(selected), args.chunk_size):
                        entries = selected[start:start+args.chunk_size]
                        buffer = np.empty((len(entries), n_samples), dtype=np.int8)
                        for i, entry in enumerate(entries):
                            buffer[i] = read_counted_dosage(reader, entry["measurement_plan"], n_samples, reader_audit)
                        for entry, counts in zip(entries, buffer):
                            for row in count_rows(entry, counts, groups, n_samples, args): cw.writerow(row)
                total["measured_nucleotides"] += len(selected)
            audit.update(chrom=c, sources=source, n_panel_samples=n_samples, n_selected_samples=len(groups[0][2]),
                         n_measurable_nucleotides=len(selected), hardcall_measurement_audit=dict(reader_audit),
                         hardcall_measurement_audit_scope="full_PSAM_before_sample_selection",
                         pgen_dimensions="checked" if selected else "not_opened_no_measurable_base", max_genotype_buffer_bytes=args.chunk_size*n_samples)
            requested = sum(len(t["requested_bases"]) for t in by_chrom[c])
            print(f"[dataset] chr={c} requested_nucleotides={requested} measurable={len(selected)} unavailable_NA={requested-len(selected)} full_PSAM_called_nucleotide_person_calls={reader_audit['called_measurements']} full_PSAM_missing_nucleotide_person_calls={reader_audit['missing_measurements']}", flush=True)
            chromosome_audits.append(audit)
    (out/"measurement_plans.json").write_text(json.dumps(dict(schema="panel-nucleotide-count-1",
        genome_build=build(args.panel_build), evidence_kind=args.evidence_kind, measurements=saved_plans), indent=2)+"\n")
    return dict(status="complete" if total["measured_nucleotides"] else "no_usable_exact_high_quality_contrast", schema=VERSION,
                counts=dict(total), target_build=build(args.target_build), panel_build=build(args.panel_build),
                target_build_contract="CLI checked against each nonmissing build/genome_build/assembly/target_build column; no relabeling or liftover",
                quality_contract=quality_proof, chromosomes=chromosome_audits, parameters=vars(args), limitations=LIMITS,
                evidence_kind=args.evidence_kind, evidence_source="user-declared; PGEN itself does not certify WGS/direct versus imputed provenance",
                genotype_semantics="PGEN hardcall 0/1/2 or missing -9; no posterior dosages; actual PSAM N; autosomal AN=2*n_called",
                memory="One PVAR pass per chromosome; compact uint64 allele offsets plus target-only metadata; bounded int8 genotype chunks, never the full variant-by-person matrix",
                inference="Prespecified nucleotide-frequency descriptions only; separate ALT measurements do not establish joint tri-allelism")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets", required=True, help="TSV: site_id chrom pos ref alts; optional requested_bases restricts counting to that REF+ALT subset")
    for name in ("target-build", "panel-build"):
        p.add_argument("--"+name, required=True, choices=("37", "38", "GRCh37", "GRCh38"))
    p.add_argument("--pgen-template", required=True, help="PGEN prefix with {chrom}, e.g. /mnt/f/gen/ukb/37/imp/chr{chrom}; .pgen/.psam appended")
    p.add_argument("--pvar-template", help="Explicit .pvar or .pvar.zst with {chrom}; required when both companions exist")
    p.add_argument("--evidence-kind", choices=("direct_gt", "imputed_hardcall"), required=True)
    p.add_argument("--quality-info", help="Exact chrom pos ref alt info table; required for imputed hardcalls")
    p.add_argument("--quality-build", choices=("37", "38", "GRCh37", "GRCh38"), help="Required for standalone INFO without the validated MFI manifest")
    p.add_argument("--min-info", type=float, default=.8)
    p.add_argument("--samples", help="Optional TSV sample_id (=unique PSAM IID), optional population; otherwise ALL uses actual PSAM roster")
    p.add_argument("--max-targets", type=int, default=5000, help="Candidate cap; 0 explicitly removes it")
    p.add_argument("--chunk-size", type=int, default=16)
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
        print(json.dumps(dict(status=result["status"], counts=result["counts"]), sort_keys=True))
        return 0
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        raise SystemExit(f"panel-targets: {exc}")


if __name__ == "__main__":
    raise SystemExit(main())
