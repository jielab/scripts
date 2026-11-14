#!/usr/bin/env python3
"""Follow up frozen different-protein hypotheses in complete signed GWAS regions.

The nominated assay keeps its true target identity. Every mapped assay of both
targets is planned before complete summaries are opened. Whole regions include
the UNION of strict source leads, not just their shared significant intersection,
plus both annotated cis windows. Shape statistics are descriptive: linked SNPs,
sample overlap, and selected regions preclude an independent replication claim.
No correlation is a colocalization probability or proof of captured identity.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import math
import os
import re
import shutil
import sqlite3
import subprocess
import tempfile
from collections import Counter, defaultdict
from decimal import Decimal
from functools import lru_cache
from pathlib import Path

import numpy as np

import assay_specificity as spec
import signed_input as si
import specificity_check as sc

SCHEMA = "phole-fingerprint-full-gwas-v1"
DATASETS = tuple(si.DATASETS)
MODES = ("FULL", "REGION_HELDOUT", "BLOCK5MB_HELDOUT", "CHROMOSOME_HELDOUT")
STANDARD = sc.STANDARD
REGION_FIELDS = ["HYPOTHESIS_ID", "REGION_ID", "BUILD", "CHR", "START", "END", "REGION_ROLES",
                 "ANCHOR_EXCLUDED", "BLOCK5MB_EXCLUDED", "CHROMOSOME_EXCLUDED", "SELECTED_SOURCE_ASSAYS"]
ASSAY_FIELDS = ["HYPOTHESIS_ID", "ROLE", "PROT", "DATASET", "ASSAY", "MAPPED", "SUMMARY_STATUS",
                "MANIFEST_STATUS", "BUILD", "ANCESTRY", "PATH", "INPUT_SCOPE", "SCOPE_EVIDENCE", "BUILD_AUTHORITY"]
LEAD_FIELDS = ["HYPOTHESIS_ID", "QUERY_ID", "REGION_ID", "BUILD", "CHR", "POS", "SNP", "PROT",
               "DATASET", "ASSAY", "REPORTED_P", "CLASS", "SELECTION_THRESHOLD", "SELECTED_SOURCE"]
DB_FIELDS = ["chrom", "pos", "snp", "ea", "nea", "beta", "se", "p", "p_text", "log10p", "log10p_text",
             "p_status", "eaf", "n", "z", "qc_reason", "flags"]
EXPORT_FIELDS = ["CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P_FLOAT", "P", "LOG10P", "LOG10P_RAW",
                 "P_STATUS", "EAF", "N", "Z", "QC_REASON", "FLAGS"]


def cell(value):
    return "." if value is None or isinstance(value, float) and not math.isfinite(value) else value


def uid(prefix, value):
    return prefix + si.digest(value)[:20]


def joined(values):
    return ";".join(sorted({str(v) for v in values if v not in (None, "", ".")})) or "."


def atomic_rows(path, rows, fields):
    si.atomic_tsv(path, ({k: cell(r.get(k)) for k in fields} for r in rows), fields)


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def protected(path, expected):
    if not expected or not Path(path).is_file() or si.file_hash(path) != expected:
        raise ValueError(f"Missing, changed, or unreceipted input: {path}")


def load_inputs(args):
    """Verify the used source files against existing cross-check receipts."""
    cross, fp = args.cross_dir, args.fingerprints_dir
    receipt_path = fp / "fingerprint-dossiers.summary.json"
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("schema") != "phole-fingerprint-dossiers-v1":
        raise ValueError("fingerprint-check requires current probe-preserving fingerprint dossiers")
    input_hashes = {"fingerprints/fingerprint-dossiers.summary.json": si.file_hash(receipt_path)}
    for name in ("cross-check.data.json", "signals.tsv", "provenance.json", "summary.json", "input-files.tsv"):
        protected(cross / name, receipt.get("input_sha256", {}).get(name))
        input_hashes[name] = si.file_hash(cross / name)
    for name in ("assay-evidence.tsv", "anchor-evidence.tsv"):
        protected(fp / name, receipt.get("output_sha256", {}).get(name))
        input_hashes["fingerprints/" + name] = si.file_hash(fp / name)
    provenance = json.loads((cross / "provenance.json").read_text())
    for name in ("coordinates-used.tsv", "prot.names.resolved.tsv"):
        protected(cross / name, provenance.get("output_sha256", {}).get(name))
        input_hashes[name] = si.file_hash(cross / name)
    payload = json.loads((cross / "cross-check.data.json").read_text())
    summary = json.loads((cross / "summary.json").read_text())
    if payload.get("generated") != provenance.get("generated") or payload.get("summary") != summary:
        raise ValueError("Cross-check data, summary, and provenance disagree")
    if payload["config"].get("strong_thresholds") != summary.get("strong_thresholds"):
        raise ValueError("Cross-check strict thresholds disagree")
    all_hypotheses = list(read_rows(fp / "assay-evidence.tsv"))
    if len({r["HYPOTHESIS_ID"] for r in all_hypotheses}) != len(all_hypotheses):
        raise ValueError("Duplicate hypothesis identity in fingerprint dossier")
    if args.pairs:
        request = list(read_rows(args.pairs))
        if not request:
            raise ValueError("--pairs must contain at least one frozen hypothesis")
        requested = {r.get("HYPOTHESIS_ID", "") for r in request}
        if "" in requested or len(requested) != len(request):
            raise ValueError("--pairs needs one unique HYPOTHESIS_ID per row from assay-evidence.tsv")
        unknown = requested - {r["HYPOTHESIS_ID"] for r in all_hypotheses}
        if unknown:
            raise ValueError("Unknown frozen HYPOTHESIS_ID: " + joined(unknown))
        hypotheses = [r for r in all_hypotheses if r["HYPOTHESIS_ID"] in requested]
        input_hashes["explicit_pairs"] = si.file_hash(args.pairs)
    else:
        hypotheses = [r for r in all_hypotheses if r["TIER"] in {"T1_distinctive_heldout", "T2_specific_full_profile"}
                      or args.selection == "all-priority" and r.get("EXACT_ANCHOR_MULTILOCUS") == "1"]
    names = {p for r in hypotheses for p in (r["PROT"], r["ALTERNATIVE_PROT"])}
    proteins = {p["prot"]: p for p in payload["proteins"] if p["prot"] in names}
    if set(proteins) != names:
        raise ValueError("Frozen target identity missing from cross-check")
    # Restrict normalization to the selected targets, but retain EVERY cohort
    # and probe for those targets. No per-hypothesis probe pooling occurs.
    signal_rows = [r for r in read_rows(cross / "signals.tsv") if r["prot"] in names]
    signals, _, _ = spec.normalize_signals(signal_rows, list(proteins.values()))
    thresholds = {d: sc.exact_p(summary["strong_thresholds"][d]) for d in DATASETS}
    if any(not Decimal(0) < p <= Decimal("5e-8") for p in thresholds.values()):
        raise ValueError("Invalid strict source threshold")
    strict = [s for s in signals if s["clean"] and s["kind"] in {"cis", "trans"}
              and s["p_decimal"] <= thresholds[s["dataset"]]]
    return hypotheses, proteins, strict, thresholds, input_hashes


def parse_windows(value):
    out = []
    for part in str(value or ".").split(";"):
        if part == ".":
            continue
        match = re.fullmatch(r"(37|38):([^:]+):(\d+)-(\d+)", part)
        if not match:
            raise ValueError("Invalid frozen exclusion interval: " + part)
        b, c, lo, hi = match.groups()
        lo, hi, c = int(lo), int(hi), si.normal_chr(c)
        if lo < 1 or hi < lo or c not in STANDARD:
            raise ValueError("Invalid frozen exclusion coordinate: " + part)
        out.append((b, c, lo, hi))
    return out


def intervals_merged(intervals):
    out = []
    for c, lo, hi in sorted(set(intervals), key=lambda x: (si.chrom_order(x[0]), x[1], x[2])):
        if out and out[-1][0] == c and lo <= out[-1][2] + 1:
            out[-1] = (c, out[-1][1], max(out[-1][2], hi))
        else:
            out.append((c, lo, hi))
    return out


def overlaps(interval, exclusions):
    b, c, lo, hi = interval
    return any(b == bb and c == cc and lo <= high and low <= hi for bb, cc, low, high in exclusions)


def build_plan(hypotheses, proteins, signals, windows, thresholds, region_flank):
    """Freeze target identities, whole windows and union leads before full input."""
    assays, regions, leads, exclusions = [], [], [], []
    for h in hypotheses:
        hid, query, alt = h["HYPOTHESIS_ID"], h["PROT"], h["ALTERNATIVE_PROT"]
        if query == alt or spec.same_coding_gene(proteins[query], proteins[alt]):
            raise ValueError("Different-protein follow-up must retain distinct coding target identities")
        for role, name in (("intended", query), ("alternative", alt)):
            protein = proteins[name]
            for d in DATASETS:
                observation = protein["data"][d]
                known = {a["alias"]: a for a in observation.get("assays", [])}
                aliases = list(observation.get("aliases", []))
                if len(aliases) != len(set(aliases)):
                    raise ValueError(f"Duplicated mapped alias for {name}/{d}")
                for alias in aliases or ["."]:
                    assays.append(dict(HYPOTHESIS_ID=hid, ROLE=role, PROT=name, DATASET=d, ASSAY=alias,
                                       MAPPED=int(alias != "."), SUMMARY_STATUS=known.get(alias, {}).get("status", "not_mapped" if alias == "." else "summary_missing")))
        expected = {(r["PROT"], r["DATASET"], r["ASSAY"]) for r in assays if r["HYPOTHESIS_ID"] == hid and r["MAPPED"]}
        if (query, h["PRIMARY"], h["PRIMARY_ASSAY"]) not in expected or (alt, h["REFERENCE"], h["REFERENCE_ASSAY"]) not in expected:
            raise ValueError("Frozen hypothesis assays disagree with their actual mapped targets")
        seeds = []
        for role, name in (("intended_cis", query), ("alternative_cis", alt)):
            w = windows.get(name)
            if not w or not w["verified"]:
                exclusions.append(dict(HYPOTHESIS_ID=hid, PROT=name, STATUS="cis_annotation_unresolved",
                                       REASON=w["reason"] if w else "annotation_missing"))
                continue
            for c, lo, hi in w["intervals"]:
                seeds.append(dict(build=w["build"], chrom=c, start=lo, end=hi, roles={role}, sources=set()))
        selected = [s for s in signals if s["prot"] in {query, alt}]
        for s in selected:
            if (s["prot"], s["dataset"], s["assay"]) not in expected:
                raise ValueError("Source lead belongs to an unmapped assay")
            seeds.append(dict(build=s["build"], chrom=s["chr"], start=max(1, s["pos"] - region_flank),
                              end=s["pos"] + region_flank, roles={"strict_lead_union"}, sources={s["dataset"] + ":" + s["assay"]}))
        merged = []
        for seed in sorted(seeds, key=lambda x: (x["build"], si.chrom_order(x["chrom"]), x["start"], x["end"])):
            if merged and (merged[-1]["build"], merged[-1]["chrom"]) == (seed["build"], seed["chrom"]) and seed["start"] <= merged[-1]["end"] + 1:
                merged[-1]["end"] = max(merged[-1]["end"], seed["end"])
                merged[-1]["roles"].update(seed["roles"])
                merged[-1]["sources"].update(seed["sources"])
            else:
                merged.append(dict(seed, roles=set(seed["roles"]), sources=set(seed["sources"])))
        anchor_windows = parse_windows(h.get("REGION_EXCLUDED"))
        blocks = parse_windows(h.get("BLOCK5MB_EXCLUDED"))
        chromosomes = {tuple(v.split(":")) for v in h.get("CHROMOSOMES_EXCLUDED", ".").split(";") if v != "."}
        for s in merged:
            interval = s["build"], s["chrom"], s["start"], s["end"]
            regions.append(dict(HYPOTHESIS_ID=hid, REGION_ID=uid("REG_", (hid, interval)), BUILD=s["build"], CHR=s["chrom"],
                                START=s["start"], END=s["end"], REGION_ROLES=joined(s["roles"]),
                                ANCHOR_EXCLUDED=int(overlaps(interval, anchor_windows)), BLOCK5MB_EXCLUDED=int(overlaps(interval, blocks)),
                                CHROMOSOME_EXCLUDED=int((s["build"], s["chrom"]) in chromosomes), SELECTED_SOURCE_ASSAYS=joined(s["sources"])))
        for s in selected:
            reg = next(r for r in regions if r["HYPOTHESIS_ID"] == hid and r["BUILD"] == s["build"]
                       and r["CHR"] == s["chr"] and r["START"] <= s["pos"] <= r["END"])
            # The identifier is a source coordinate/rsID query, not a guessed
            # allele. Actual EA/NEA are obtained from a full source record.
            qid = uid("VAR_", (s["build"], s["chr"], s["pos"], s["snp"]))
            leads.append(dict(HYPOTHESIS_ID=hid, QUERY_ID=qid, REGION_ID=reg["REGION_ID"], BUILD=s["build"], CHR=s["chr"], POS=s["pos"],
                              SNP=s["snp"], PROT=s["prot"], DATASET=s["dataset"], ASSAY=s["assay"], REPORTED_P=s["p"],
                              CLASS=s["kind"], SELECTION_THRESHOLD=str(thresholds[s["dataset"]]), SELECTED_SOURCE=1))
    return assays, regions, leads, exclusions


def exact_cis_plan(hypotheses, windows):
    """Keep genuine gene boundaries separate from merged display regions."""
    rows = []
    for h in hypotheses:
        for role, target in (("intended_cis", h["PROT"]), ("alternative_cis", h["ALTERNATIVE_PROT"])):
            w = windows.get(target, {})
            intervals = w.get("intervals", []) if w.get("verified") else []
            for i, (chrom, lo, hi) in enumerate(intervals or [(".", None, None)], 1):
                rows.append(dict(HYPOTHESIS_ID=h["HYPOTHESIS_ID"], WINDOW_ID=uid("CIS_", (h["HYPOTHESIS_ID"], target, w.get("build"), chrom, lo, hi)),
                                 WINDOW_PROT=target, WINDOW_ROLE=role, COMPONENT_INDEX=i if intervals else 0, N_COMPONENTS=len(intervals),
                                 ANNOTATION_VERIFIED=int(bool(w.get("verified") and intervals)), BUILD=w.get("build") or ".",
                                 CHR=chrom, START=lo, END=hi, REASON=w.get("reason", "annotation_missing")))
    return rows


def exact_cis_checks(plan, assays, dbs, receipts, thresholds):
    rows = []
    for w in plan:
        for a in (r for r in assays if r["HYPOTHESIS_ID"] == w["HYPOTHESIS_ID"]):
            region = dict(REGION_ID=w["WINDOW_ID"], BUILD=w["BUILD"], CHR=w["CHR"], START=w["START"], END=w["END"], REGION_ROLES=w["WINDOW_ROLE"])
            inspection_assay = a if w["ANNOTATION_VERIFIED"] else dict(a, MANIFEST_STATUS="cis_annotation_unresolved")
            key = a["DATASET"], a["ASSAY"]
            source_rows = fetch_region(dbs.get(key), region) if region_status(inspection_assay, region) == "ready" else []
            qc = region_qc(w["HYPOTHESIS_ID"], region, inspection_assay, source_rows, receipts.get(key), thresholds[a["DATASET"]])
            rows.append(dict(qc, WINDOW_ID=w["WINDOW_ID"], WINDOW_PROT=w["WINDOW_PROT"], WINDOW_ROLE=w["WINDOW_ROLE"],
                             COMPONENT_INDEX=w["COMPONENT_INDEX"], N_COMPONENTS=w["N_COMPONENTS"], ANNOTATION_VERIFIED=w["ANNOTATION_VERIFIED"],
                             IS_ASSAY_OWN_CIS=int(a["PROT"] == w["WINDOW_PROT"]), BOUNDARY_SCOPE="exact_annotated_gene_body_plus_cis_flank_not_merged_lead_region"))
    groups = defaultdict(list)
    for r in rows:
        groups[r["HYPOTHESIS_ID"], r["WINDOW_PROT"], r["PROT"], r["DATASET"], r["ASSAY"]].append(r)
    for rr in groups.values():
        complete = bool(rr[0]["N_COMPONENTS"] and len(rr) == rr[0]["N_COMPONENTS"] and all(r["DECLARED_REGION_READ_COMPLETE"] for r in rr))
        any_gws = any(r["HAS_GWS"] == 1 for r in rr)
        has_gws = 1 if any_gws else 0 if complete and all(r["HAS_GWS"] == 0 for r in rr) else None
        for r in rr:
            r.update(WHOLE_GENE_WINDOW_READ_COMPLETE=int(complete), WHOLE_GENE_WINDOW_HAS_GWS=has_gws,
                     WHOLE_GENE_WINDOW_NO_GWS_OBSERVED=int(has_gws == 0) if complete and has_gws is not None else None)
    return rows


def build_authority(row, manifest_path):
    """Reject legacy cross-check inference as sole full-file build evidence.

    An explicit user manifest without generated evidence fields is a build
    declaration. When it does name a readable receipt, its declared scope and
    fingerprint are retained; a dataset-level summary inference is not silently
    upgraded to knowledge about the complete file.
    """
    evidence = row.get("BUILD_EVIDENCE", "")
    if evidence in {"", "."}:
        return "explicit_manifest_build_declaration", []
    sources, inferred, authoritative, conflict = [], False, False, False
    full_path = Path(row["PATH"]).resolve()

    def agrees(raw):
        nonlocal conflict
        try:
            build = si.normal_build(raw)
        except ValueError:
            return False
        if build != row["BUILD"]:
            conflict = True
            return False
        return True

    def matching_path(raw, base):
        p = Path(os.path.expandvars(os.path.expanduser(raw)))
        return (base / p if not p.is_absolute() else p).resolve() == full_path
    for item in evidence.split(";"):
        if item.strip() == "explicit --dataset-builds":
            authoritative = True
            continue
        source, _, selector = item.partition("#")
        path = Path(source)
        if not path.is_absolute():
            path = Path(manifest_path).resolve().parent / path
        if not path.is_file():
            continue
        sources.append(dict(path=str(path.resolve()), sha256=si.file_hash(path)))
        if selector.startswith(("DATASET=", "ASSAY=")):
            column, value = selector.split("=", 1)
            selected = [{k.upper(): v for k, v in r.items()} for r in read_rows(path)]
            selected = [r for r in selected if ("ukb" if r.get(column) == "prot" else r.get(column)) == value] if column == "DATASET" else [r for r in selected if r.get(column) == value]
            for r in selected:
                if r.get("PATH") not in (None, "", ".") and not matching_path(r["PATH"], path.parent):
                    continue
                if not agrees(r.get("BUILD") or r.get("INFERRED_BUILD") or r.get("BEST_BUILD")):
                    continue
                if "inference" in r.get("SCOPE", "").lower() or "does not certify" in r.get("SCOPE", "").lower():
                    inferred = True
                else:
                    authoritative = True
        else:
            # Recheck both freshness AND content when consuming an older
            # generated manifest. A fresh receipt may now describe a new build.
            fresh = full_path.is_file() and path.stat().st_mtime_ns >= full_path.stat().st_mtime_ns
            if path.suffix == ".grch":
                tokens = path.read_text().strip().split()
                valid = len(tokens) == 1 and agrees(tokens[0])
                if fresh and valid:
                    authoritative = True
            elif path.name.endswith(".sig.meta.tsv"):
                meta = {r.get("key", "").lower(): r.get("value", "") for r in read_rows(path)}
                if meta.get("association_scope") == "full_gwas" and meta.get("association_input") and Path(meta["association_input"]).is_absolute() and matching_path(meta["association_input"], path.parent):
                    valid = agrees(meta.get("build"))
                    if fresh and valid:
                        authoritative = True
    if conflict:
        return "full_file_build_receipt_conflict", sources
    if inferred and not authoritative:
        return "summary_build_inference_does_not_certify_full_file", sources
    return "recorded_manifest_build_evidence" if authoritative else "full_file_build_evidence_stale_or_unavailable", sources


def attach_manifest(plan, path):
    """Keep mapped-but-missing entries; reject deliberate target relabelling."""
    expected = {}
    for row in plan:
        if not row["MAPPED"]:
            continue
        key = row["DATASET"], row["ASSAY"]
        if key in expected and expected[key] != row["PROT"]:
            raise ValueError("A mapped assay has two canonical identities: " + str(key))
        expected[key] = row["PROT"]
    manifests = {}
    seen = set()
    path_owners, physical_owners = {}, {}
    for raw in read_rows(path) if path else ():
        if not {"DATASET", "ASSAY", "PROT", "PATH", "BUILD", "INPUT_SCOPE"}.issubset(raw):
            raise ValueError("Manifest requires DATASET ASSAY PROT PATH BUILD INPUT_SCOPE")
        r = {k: str(v).strip() for k, v in raw.items() if k is not None}
        d, alias = sc.dataset(r["DATASET"]), si.safe_assay(r["ASSAY"])
        key = d, alias
        if key in seen:
            raise ValueError("Duplicate dataset/assay in manifest: " + str(key))
        seen.add(key)
        if key not in expected:
            continue
        if r["PROT"] != expected[key]:
            raise ValueError(f"Manifest target identity mismatch for {d}/{alias}: {r['PROT']} != {expected[key]}")
        if r["INPUT_SCOPE"].lower() != "full":
            raise ValueError(f"{d}/{alias}: declared INPUT_SCOPE=full is required")
        p = Path(os.path.expandvars(os.path.expanduser(r["PATH"])))
        if not p.is_absolute():
            p = Path(path).resolve().parent / p
        p = p.resolve()
        if p in path_owners and path_owners[p] != key:
            raise ValueError(f"Same full GWAS file reused for distinct assay identities: {path_owners[p]} and {key}: {p}")
        path_owners[p] = key
        if p.is_file():
            stat = p.stat()
            # st_ino=0 is not a usable physical identity on some mounted file
            # systems. Nonzero device/inode identifies symlinks or hardlinks
            # without a prohibitively expensive whole-GWAS content checksum.
            physical = stat.st_dev, stat.st_ino
            if stat.st_ino and physical in physical_owners and physical_owners[physical] != key:
                raise ValueError(f"Same full GWAS file reused for distinct assay identities (physical file): {physical_owners[physical]} and {key}: {p}")
            if stat.st_ino:
                physical_owners[physical] = key
        if p.name.endswith((".sig.txt", ".genes.out", ".genes.raw")) or re.search(r"(?:^|[._-])(hm3|cis|thin|lead|sig|magma)(?:[._-]|$)", p.name, re.I):
            raise ValueError("Filtered/selected input is not a full GWAS: " + str(p))
        r.update(DATASET=d, ASSAY=alias, PATH=str(p), ANCESTRY=r.get("ANCESTRY", "UNKNOWN").upper() or "UNKNOWN",
                 EFFECT_SCALE=r.get("EFFECT_SCALE", "unknown"), SYNTHETIC=r.get("SYNTHETIC", "FALSE").upper())
        if r["SYNTHETIC"] not in {"TRUE", "FALSE"}:
            raise ValueError("SYNTHETIC must be TRUE or FALSE")
        try:
            r["BUILD"] = si.normal_build(r["BUILD"])
            r["STATUS"] = "ready" if p.is_file() else "full_file_missing"
        except ValueError:
            r["BUILD"], r["STATUS"] = ".", "manifest_build_unresolved"
        r["BUILD_AUTHORITY"], r["BUILD_RECEIPT_FINGERPRINTS"] = build_authority(r, path)
        if r["STATUS"] == "ready" and r["BUILD_AUTHORITY"] == "summary_build_inference_does_not_certify_full_file":
            r["STATUS"] = "full_file_build_unverified_summary_inference"
        elif r["STATUS"] == "ready" and r["BUILD_AUTHORITY"] == "full_file_build_evidence_stale_or_unavailable":
            r["STATUS"] = "full_file_build_evidence_stale_or_unavailable"
        elif r["STATUS"] == "ready" and r["BUILD_AUTHORITY"] == "full_file_build_receipt_conflict":
            r["STATUS"] = "full_file_build_receipt_conflict"
        if r.get("INDEX") not in (None, "", "."):
            index = Path(r["INDEX"])
            r["INDEX"] = str((Path(path).resolve().parent / index if not index.is_absolute() else index).resolve())
        else:
            r.pop("INDEX", None)
        manifests[key] = r
    synthetic = {r["SYNTHETIC"] for r in manifests.values()}
    if len(synthetic) > 1:
        raise ValueError("Synthetic and real full inputs cannot be mixed")
    for row in plan:
        m = manifests.get((row["DATASET"], row["ASSAY"]))
        status = "not_mapped" if not row["MAPPED"] else "plan_only_not_inspected" if path is None else "manifest_assay_missing"
        row.update(MANIFEST_STATUS=m["STATUS"] if m else status, BUILD=m["BUILD"] if m else ".",
                   ANCESTRY=m["ANCESTRY"] if m else "UNKNOWN", PATH=m["PATH"] if m else ".",
                   INPUT_SCOPE=m["INPUT_SCOPE"] if m else ".", SCOPE_EVIDENCE=m.get("SCOPE_EVIDENCE", "explicit_manifest_declaration") if m else ".",
                   BUILD_AUTHORITY=m["BUILD_AUTHORITY"] if m else ".")
    return manifests, next(iter(synthetic), None) == "TRUE" if synthetic else None


def full_record(fields, header, default_n=None):
    raw = dict(zip(header, fields))
    c, pos = si.normal_chr(raw["CHR"]), int(raw["POS"])
    p_text, reason, flags = raw.get("P", "."), ".", []
    try:
        precise = sc.exact_p(p_text)
        p = float(precise)
        p_status = "source_zero_underflow" if precise == 0 else "nonzero_source_float_underflow" if p == 0 else "positive_source_P"
    except ValueError:
        precise, p, p_status = None, None, "invalid_source_P"
    log_text, logp = raw.get("LOG10P", "."), None
    if log_text.strip().upper() not in si.MISSING:
        try:
            logp = si.number(log_text)
            if logp < 0:
                raise ValueError()
            if precise and abs(float(-precise.log10()) - logp) > 0.05:
                flags.append("source_LOG10P_P_inconsistent")
                logp = None
        except (ValueError, OverflowError):
            flags.append("invalid_source_LOG10P")
            logp = None
    if logp is None and precise and precise > 0:
        logp = float(-precise.log10())
    if precise == 0:
        # P may be rounded, clipped, or have underflowed upstream. A supplied
        # finite -log10(P) is more informative than literal zero; do not invent
        # a 1e-300 bound, nor call a moderate valid LOG10P corrupt merely because
        # P was printed as zero.
        p_status = "source_zero_LOG10P_available" if logp is not None else "source_zero_precision_unknown"
        flags.append("literal_source_zero_not_exact_probability")
    try:
        r = si.parse_record(fields, header, default_n)
        if precise is None:
            raise ValueError("invalid_P")
        vals = [r["beta"], r["se"], r["eaf"], r["n"], r["z"]]
        vals = [v if math.isfinite(v) else None for v in vals]
        beta, se, eaf, n, z = vals
    except ValueError as error:
        reason = str(error)
        beta = se = eaf = n = z = None
    return (c, pos, raw.get("SNP", "."), raw.get("EA", ".").upper(), raw.get("NEA", ".").upper(),
            beta, se, p, p_text, logp, log_text, p_status, eaf, n, z, reason, joined(flags))


def in_intervals(c, p, bounds):
    if c not in bounds:
        return False
    starts, intervals = bounds[c]
    i = bisect.bisect_right(starts, p) - 1
    return i >= 0 and p <= intervals[i][1]


def extract_regions_cached(assay, intervals, args, implementation):
    """Scan each complete file once; spill regional rows to indexed SQLite."""
    intervals = intervals_merged(intervals)
    fingerprint = si.summary_identity(assay)
    signature = dict(schema=SCHEMA, input=fingerprint, assay=assay, intervals=intervals,
                     implementation=implementation, plan_sha256=getattr(args, "plan_sha256", None),
                     manifest_sha256=getattr(args, "manifest_sha256", None))
    key = si.digest(signature)
    stem = assay["DATASET"] + "-" + si.digest(assay["ASSAY"])[:20]
    state = args.out / "state"
    export = args.out / "extracted" / (stem + ".tsv.gz")
    state.mkdir(parents=True, exist_ok=True)
    export.parent.mkdir(parents=True, exist_ok=True)
    db_path, receipt_path = state / (stem + ".sqlite3"), state / (stem + ".json")
    if receipt_path.is_file():
        saved = json.loads(receipt_path.read_text())
        if saved.get("checkpoint_key") == key and db_path.is_file() and export.is_file():
            if si.file_hash(db_path) != saved.get("database_sha256") or si.file_hash(export) != saved.get("extracted_sha256"):
                raise ValueError("Changed regional checkpoint; remove its state files before rerunning: " + str(db_path))
            return db_path, dict(saved, checkpoint_reused=True)
    path = Path(assay["PATH"])
    header = si.read_header(path)
    ci, pi = header.index("CHR"), header.index("POS")
    grouped = defaultdict(list)
    for c, lo, hi in intervals:
        grouped[c].append((lo, hi))
    bounds = {c: ([x[0] for x in values], values) for c, values in grouped.items()}
    stats = Counter()
    fd, temp_name = tempfile.mkstemp(prefix=".fingerprint-", suffix=".sqlite3", dir=state)
    os.close(fd)
    temporary_db = Path(temp_name)
    temporary_export = export.with_name("." + export.name + ".tmp")
    con = sqlite3.connect(temporary_db)
    con.execute("CREATE TABLE variants (chrom TEXT, pos INTEGER, snp TEXT, ea TEXT, nea TEXT, beta REAL, se REAL, p REAL, p_text TEXT, log10p REAL, log10p_text TEXT, p_status TEXT, eaf REAL, n REAL, z REAL, qc_reason TEXT, flags TEXT)")
    buffer = []

    def consume(line):
        if not line.strip() or line.startswith("#"):
            return
        fields = line.split()
        if len(fields) > 1 and fields[0].upper() == header[0] and fields[1].upper() == header[1]:
            return
        stats["rows_read"] += 1
        try:
            if len(fields) != len(header):
                raise ValueError()
            c, pos = si.normal_chr(fields[ci]), int(fields[pi])
            if c not in STANDARD or pos < 1:
                raise ValueError()
        except (ValueError, IndexError):
            stats["unlocatable_rows"] += 1
            return
        if not in_intervals(c, pos, bounds):
            return
        row = full_record(fields, header, assay.get("N"))
        stats["regional_rows"] += 1
        stats["regional_valid_effect_rows" if row[-2] == "." else "regional_failed_QC_rows"] += 1
        if row[7] is not None:
            stats["regional_valid_P_rows"] += 1
        if row[11] in {"source_zero_LOG10P_available", "source_zero_precision_unknown", "nonzero_source_float_underflow"}:
            stats[row[11]] += 1
        buffer.append(row)
        if len(buffer) >= 4096:
            con.executemany("INSERT INTO variants VALUES (" + ",".join("?" for _ in DB_FIELDS) + ")", buffer)
            buffer.clear()

    index_paths = [Path(assay["INDEX"])] if assay.get("INDEX") else [Path(str(path) + ".tbi"), Path(str(path) + ".csi")]
    index = next((p for p in index_paths if p.is_file()), None)
    executable = shutil.which(args.tabix)
    method = "explicit_full_file_scan"
    try:
        if assay.get("INDEX") and index is None:
            raise ValueError("Explicit INDEX is missing: " + assay["INDEX"])
        if index is not None and executable is not None:
            indexed_path = str(path) if index in (Path(str(path) + ".tbi"), Path(str(path) + ".csi")) else str(path) + "##idx##" + str(index)
            listed = subprocess.run([executable, "-l", indexed_path], check=True, text=True, capture_output=True).stdout.splitlines()
            contigs = defaultdict(list)
            for c in listed:
                contigs[si.normal_chr(c)].append(c)
            if any(len(v) != 1 for v in contigs.values()):
                raise ValueError("Ambiguous normalized contigs in tabix index")
            with tempfile.TemporaryDirectory(prefix="phole-fingerprint-") as temporary:
                bed, err = Path(temporary) / "regions.tsv", Path(temporary) / "stderr"
                with bed.open("w") as handle:
                    for c, lo, hi in intervals:
                        if c in contigs:
                            handle.write(f"{contigs[c][0]}\t{lo}\t{hi}\n")
                        else:
                            stats["requested_intervals_absent_contig"] += 1
                if bed.stat().st_size:
                    with err.open("w") as errors:
                        process = subprocess.Popen([executable, "-R", str(bed), indexed_path], stdout=subprocess.PIPE, stderr=errors, text=True)
                        try:
                            for line in process.stdout:
                                consume(line)
                        finally:
                            process.stdout.close()
                            rc = process.wait()
                    if rc:
                        raise RuntimeError("tabix failed; scan was not substituted: " + err.read_text()[-2000:])
            method = "tabix_complete_regions"
        else:
            if not args.allow_scan:
                raise ValueError(f"{path}: usable tabix index/executable required or explicit --allow-scan")
            with si.open_text(path) as handle:
                for line in handle:
                    consume(line)
        if buffer:
            con.executemany("INSERT INTO variants VALUES (" + ",".join("?" for _ in DB_FIELDS) + ")", buffer)
        con.execute("CREATE INDEX loci ON variants (chrom,pos)")
        con.commit()
        with gzip.open(temporary_export, "wt", newline="", encoding="utf-8") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(EXPORT_FIELDS)
            for row in con.execute("SELECT * FROM variants ORDER BY chrom,pos,snp,ea,nea"):
                writer.writerow([cell(v) for v in row])
        con.close()
        if si.summary_identity(assay) != fingerprint:
            raise ValueError("Complete GWAS changed during regional extraction: " + str(path))
        os.replace(temporary_db, db_path)
        os.replace(temporary_export, export)
        receipt = dict(DATASET=assay["DATASET"], ASSAY=assay["ASSAY"], PROT=assay["PROT"], BUILD=assay["BUILD"],
                       SYNTHETIC=assay["SYNTHETIC"], input_fingerprint=fingerprint, intervals=intervals, method=method,
                       statistics=dict(stats), header=header, index=str(index) if index else None, checkpoint_key=key,
                       checkpoint_reused=False, database=str(db_path.relative_to(args.out)), database_sha256=si.file_hash(db_path),
                       extracted=str(export.relative_to(args.out)), extracted_sha256=si.file_hash(export))
        si.atomic_json(receipt_path, receipt)
        return db_path, receipt
    finally:
        con.close()
        temporary_db.unlink(missing_ok=True)
        temporary_export.unlink(missing_ok=True)


def fetch_region(db, region):
    if db is None:
        return []
    with sqlite3.connect(db) as con:
        raw = con.execute("SELECT * FROM variants WHERE chrom=? AND pos BETWEEN ? AND ? ORDER BY pos,snp,ea,nea",
                          (region["CHR"], region["START"], region["END"])).fetchall()
    records = []
    for values in raw:
        r = dict(zip(DB_FIELDS, values))
        r["chr"] = r.pop("chrom")
        r["eaf"], r["n"] = (math.nan if r[k] is None else r[k] for k in ("eaf", "n"))
        records.append(r)
    return records


def strict_unique(records):
    row, reason = sc.strict_unique(records)
    if row is None:
        return row, reason
    for other in records[1:]:
        for field in ("beta", "se"):
            if not math.isclose(row[field], other[field], rel_tol=1e-7, abs_tol=0):
                return None, "duplicate_precise_effect_conflict"
    # Underflowed P=0 records with different reported -log10(P) are not
    # interchangeable literal repeats.
    observed = [r.get("log10p") for r in records if r.get("log10p") is not None]
    if observed and max(observed) - min(observed) > 0.05:
        return None, "duplicate_LOG10P_conflict"
    return row, reason


def allele_group(records):
    grouped = defaultdict(list)
    for r in records:
        if r["qc_reason"] == ".":
            grouped[(r["pos"], sc.allele_pair_key(r["ea"], r["nea"]))].append(r)
    result, rejected = {}, Counter()
    for key, rr in grouped.items():
        r, why = strict_unique(rr)
        if r is None:
            rejected[why] += 1
        else:
            result[key] = r
    return result, rejected


def anchor_from(row, assay):
    return dict(CHR=row["chr"], POS=row["pos"], EA=row["ea"], NEA=row["nea"], EAF=row["eaf"], ANCESTRY=assay.get("ANCESTRY", "UNKNOWN"))


def correlation(x, y):
    if len(x) < 3:
        return None
    a, b = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)) or np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    value = float(np.corrcoef(a, b)[0, 1])
    return value if math.isfinite(value) else None


def ranks(values):
    values = np.asarray(values)
    order = np.argsort(values, kind="mergesort")
    out = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        out[order[i:j]] = (i + j - 1) / 2
        i = j
    return out


def shape_statistics(left, right):
    return dict(N_ALIGNED=len(left), PEARSON_Z=correlation(left, right), SPEARMAN_Z=correlation(ranks(left), ranks(right)),
                SAME_SIGN_FRACTION=sum(x * y > 0 for x, y in zip(left, right)) / len(left) if left else None)


def assay_key(row):
    return row["PROT"], row["DATASET"], row["ASSAY"]


def region_status(assay, region):
    if assay["MANIFEST_STATUS"] != "ready":
        return assay["MANIFEST_STATUS"]
    if assay["BUILD"] != region["BUILD"]:
        return "build_mismatch"
    return "ready"


def comparison_specs(h, assays):
    """Target and reciprocal controls preserve canonical identities and probes."""
    rows = [a for a in assays if a["HYPOTHESIS_ID"] == h["HYPOTHESIS_ID"]]
    primary = next(a for a in rows if a["PROT"] == h["PROT"] and a["DATASET"] == h["PRIMARY"] and a["ASSAY"] == h["PRIMARY_ASSAY"])
    alternative_primary = [a for a in rows if a["ROLE"] == "alternative" and a["DATASET"] == h["PRIMARY"]]
    comparisons = {}
    def add(left, right, role):
        if assay_key(left) == assay_key(right):
            return
        key = assay_key(left), assay_key(right)
        if key not in comparisons:
            comparisons[key] = (left, right, set())
        comparisons[key][2].add(role)
    for a in rows:
        add(primary, a, "intended_target_control" if a["ROLE"] == "intended" else "alternative_target_comparison")
        for alt in alternative_primary:
            add(alt, a, "alternative_target_control" if a["ROLE"] == "alternative" else "reciprocal_target_control")
    # In particular expose CKB Q-versus-ALT directly. Neither may be the
    # nominating primary assay, and an apparent raw-sign discrepancy must be
    # resolved using their own allele columns within that cohort.
    for d in DATASETS:
        for intended in (a for a in rows if a["DATASET"] == d and a["ROLE"] == "intended"):
            for alt in (a for a in rows if a["DATASET"] == d and a["ROLE"] == "alternative"):
                add(intended, alt, "within_cohort_target_contrast")
    triples = []
    for d in DATASETS:
        if d == h["PRIMARY"]:
            continue
        for intended in (a for a in rows if a["DATASET"] == d and a["ROLE"] == "intended"):
            for alt in (a for a in rows if a["DATASET"] == d and a["ROLE"] == "alternative"):
                triples.append((primary, intended, alt))
    return list(comparisons.values()), triples


def compare_region(left_rows, right_rows, left_assay, right_assay):
    left, lbad = allele_group(left_rows)
    right, rbad = allele_group(right_rows)
    common = set(left) & set(right)
    x, y, aligned, rejected = [], [], [], Counter(lbad + rbad)
    for key in sorted(common):
        l, r = left[key], right[key]
        # rsIDs can be missing/renamed; coordinate plus explicit alleles is the
        # variant key. Discordant extant rsIDs remain an explicit warning.
        anchor = anchor_from(l, left_assay)
        beta, eaf, reason = si.align_record(r, anchor, right_assay["ANCESTRY"])
        if beta is None:
            rejected[reason] += 1
            continue
        lz, rz = l["beta"] / l["se"], beta / r["se"]
        x.append(lz); y.append(rz)
        aligned.append(dict(CHR=l["chr"], POS=l["pos"], EA=l["ea"], NEA=l["nea"], LEFT_SNP=l["snp"], RIGHT_SNP=r["snp"],
                            LEFT_Z=lz, RIGHT_Z=rz, LEFT_P=l["p_text"], RIGHT_P=r["p_text"], ALIGNMENT=reason,
                            IDENTIFIER_CONFLICT=int(bool(si.snp_tokens(l["snp"]) and si.snp_tokens(r["snp"]) and not si.snp_tokens(l["snp"]) & si.snp_tokens(r["snp"])))))
    return dict(N_SHARED_COORDINATES=len({k[0] for k in common}), N_SHARED_ALLELE_VARIANTS=len(common),
                N_LEFT_VALID_VARIANTS=len(left), N_RIGHT_VALID_VARIANTS=len(right), N_UNRESOLVED=sum(rejected.values()),
                REJECTIONS=json.dumps(dict(rejected), sort_keys=True), **shape_statistics(x, y)), aligned


def triple_region(primary_rows, intended_rows, alternative_rows, primary, intended, alternative):
    p, pb = allele_group(primary_rows)
    q, qb = allele_group(intended_rows)
    a, ab = allele_group(alternative_rows)
    common = set(p) & set(q) & set(a)
    zp, zq, za, rejected = [], [], [], Counter(pb + qb + ab)
    for key in sorted(common):
        pr, qr, ar = p[key], q[key], a[key]
        anchor = anchor_from(pr, primary)
        bq, _, rq = si.align_record(qr, anchor, intended["ANCESTRY"])
        ba, _, ra = si.align_record(ar, anchor, alternative["ANCESTRY"])
        if bq is None or ba is None:
            rejected["intended:" + rq if bq is None else "alternative:" + ra] += 1
            continue
        zp.append(pr["beta"] / pr["se"]); zq.append(bq / qr["se"]); za.append(ba / ar["se"])
    own, alt = shape_statistics(zp, zq), shape_statistics(zp, za)
    return dict(N_COMMON_ALLELE_VARIANTS=len(common), N_COMMON_ALIGNED=len(zp),
                INTENDED_PEARSON_Z=own["PEARSON_Z"], ALTERNATIVE_PEARSON_Z=alt["PEARSON_Z"],
                DELTA_PEARSON_Z=alt["PEARSON_Z"] - own["PEARSON_Z"] if alt["PEARSON_Z"] is not None and own["PEARSON_Z"] is not None else None,
                INTENDED_SPEARMAN_Z=own["SPEARMAN_Z"], ALTERNATIVE_SPEARMAN_Z=alt["SPEARMAN_Z"],
                DELTA_SPEARMAN_Z=alt["SPEARMAN_Z"] - own["SPEARMAN_Z"] if alt["SPEARMAN_Z"] is not None and own["SPEARMAN_Z"] is not None else None,
                INTENDED_SAME_SIGN=own["SAME_SIGN_FRACTION"], ALTERNATIVE_SAME_SIGN=alt["SAME_SIGN_FRACTION"],
                REJECTIONS=json.dumps(dict(rejected), sort_keys=True))


def association_p(record):
    """Finite observed association probability, with literal-zero ambiguity."""
    p = sc.exact_p(record["p_text"])
    if p > 0:
        return p
    if record.get("log10p") is not None:
        return Decimal(10) ** (-Decimal(str(record["log10p"])))
    return None


def region_qc(hid, region, assay, rows, receipt, strict_threshold):
    valid = [r for r in rows if r["qc_reason"] == "."]
    readable = [r for r in rows if r["p"] is not None]
    # Literal source zeros are reported as underflow, with no fabricated exact
    # strength or ordering among such rows.
    best = min(readable, key=lambda r: (sc.exact_p(r["p_text"]), r["pos"], r["snp"])) if readable else None
    accepted, duplicates = allele_group(valid)
    accepted_rows = list(accepted.values())
    accepted_p = [(r, association_p(r)) for r in accepted_rows]
    unknown_precision = sum(p is None for r, p in accepted_p)
    resolved_p = [(r, p) for r, p in accepted_p if p is not None]
    accepted_best = min(resolved_p, key=lambda x: (x[1], x[0]["pos"], x[0]["snp"]))[0] if resolved_p else None
    status = region_status(assay, region)
    if status == "ready":
        status = "observed" if valid else "region_rows_all_failed_QC" if rows else "region_without_variants"
    unlocatable = receipt.get("statistics", {}).get("unlocatable_rows", 0) if receipt else 0
    complete = status == "observed" and len(valid) == len(rows) and not unlocatable and not duplicates
    def significant(threshold):
        if any(p <= threshold for r, p in resolved_p):
            return True
        return False if resolved_p and not unknown_precision and not duplicates else None
    gws, strict = significant(Decimal("5e-8")), significant(strict_threshold)
    return dict(HYPOTHESIS_ID=hid, REGION_ID=region["REGION_ID"], PROT=assay["PROT"], DATASET=assay["DATASET"], ASSAY=assay["ASSAY"],
                ROLE=assay["ROLE"], BUILD=region["BUILD"], CHR=region["CHR"], START=region["START"], END=region["END"],
                REGION_ROLES=region["REGION_ROLES"], STATUS=status, N_ROWS=len(rows), N_VALID_EFFECT_ROWS=len(valid),
                N_VALID_P_ROWS=len(readable), N_POSITIONS=len({r["pos"] for r in rows}), N_QC_FAILED_ROWS=len(rows) - len(valid),
                N_DUPLICATE_CONFLICTS=sum(duplicates.values()), UNLOCATABLE_INPUT_ROWS=unlocatable,
                N_RESOLVED_ALLELE_VARIANTS=len(accepted_rows), N_P_PRECISION_UNRESOLVED=unknown_precision,
                MIN_P=best["p_text"] if best else None, MIN_P_SNP=best["snp"] if best else None,
                MIN_P_POS=best["pos"] if best else None, MIN_P_STATUS=best["p_status"] if best else None,
                MIN_P_LOG10P=best["log10p"] if best else None, MIN_P_QC_REASON=best["qc_reason"] if best else None,
                MIN_RESOLVED_P=accepted_best["p_text"] if accepted_best else None,
                MIN_RESOLVED_SNP=accepted_best["snp"] if accepted_best else None,
                MIN_RESOLVED_LOG10P=accepted_best["log10p"] if accepted_best else None,
                HAS_GWS=None if gws is None else int(gws),
                HAS_STRICT=None if strict is None else int(strict), STRICT_THRESHOLD=str(strict_threshold),
                DECLARED_REGION_READ_COMPLETE=int(complete), NO_GWS_OBSERVED=int(not gws) if complete and gws is not None else None,
                INPUT_SCOPE=assay["INPUT_SCOPE"],
                INTERPRETATION="Association observed among valid supplied variants" if gws else
                "No genome-wide significant variant among fully read declared-region rows; this is not proof of no effect or complete population variant coverage" if complete and gws is False else
                "Absence unresolved: mapping, file, build, coverage or row QC unavailable")


def plot_profile(assay, rows, status, maximum=500):
    valid = [r for r in rows if r["qc_reason"] == "."]
    if len(valid) > maximum:
        indices = set(np.linspace(0, len(valid) - 1, max(1, maximum - 1), dtype=int))
        indices.add(max(range(len(valid)), key=lambda i: abs(valid[i]["z"])))
        displayed = [valid[i] for i in sorted(indices)]
    else:
        displayed = valid
    fields = ("chr", "pos", "snp", "ea", "nea", "z", "p_text", "log10p", "log10p_text", "p_status")
    return dict(PROT=assay["PROT"], DATASET=assay["DATASET"], ASSAY=assay["ASSAY"], ROLE=assay["ROLE"], STATUS=status,
                N_TOTAL_VALID_VARIANTS=len(valid), N_DISPLAY_VARIANTS=len(displayed), DISPLAY_SAMPLED=len(displayed) < len(valid),
                points=[{k: r[k] for k in fields} for r in displayed])


def source_record(records, snp):
    good = [r for r in records if r["qc_reason"] == "."]
    selected = [r for r in good if r["snp"] == snp or si.snp_tokens(r["snp"]) & si.snp_tokens(snp)]
    if selected:
        return strict_unique(selected)
    if len(good) == 1:
        if si.snp_tokens(snp) and si.snp_tokens(good[0]["snp"]):
            return None, "source_variant_identifier_conflict"
        return good[0], "source_coordinate_ID_differs"
    return strict_unique(good)


def compare_reported_p(reported, current, threshold):
    if current is None:
        return "source_unresolved"
    old, new = sc.exact_p(reported), sc.exact_p(current["p_text"])
    effective = association_p(current)
    if effective is None:
        return "source_full_P_zero_precision_unresolved"
    if effective > sc.exact_p(threshold):
        return "source_full_P_fails_frozen_selection_threshold"
    if old == 0 or new == 0:
        return "source_underflow_not_exactly_comparable"
    if abs(old - new) > max(old, new) * Decimal("0.01"):
        return "source_summary_P_mismatch"
    return "source_summary_P_agrees_within_1pct"


def lead_family(assays, leads, alpha):
    sources = {(r["BUILD"], r["CHR"], r["POS"], r["DATASET"], r["ASSAY"]) for r in leads}
    family = {}
    for q in leads:
        for a in assays:
            if a["HYPOTHESIS_ID"] != q["HYPOTHESIS_ID"] or not a["MAPPED"]:
                continue
            key = q["QUERY_ID"], a["PROT"], a["DATASET"], a["ASSAY"]
            if (q["BUILD"], q["CHR"], q["POS"], a["DATASET"], a["ASSAY"]) in sources:
                continue
            family[key] = dict(TEST_ID=uid("TEST_", key), QUERY_ID=q["QUERY_ID"], PROT=a["PROT"], DATASET=a["DATASET"],
                               ASSAY=a["ASSAY"], BUILD=q["BUILD"], CHR=q["CHR"], POS=q["POS"], SNP=q["SNP"],
                               FAMILY_ROLE="unselected_at_this_coordinate_exploratory_followup")
    threshold = Decimal(str(alpha)) / len(family) if family else None
    for r in family.values():
        r.update(FAMILY_SIZE=len(family), FOLLOWUP_THRESHOLD=str(threshold) if threshold else ".")
    return family, sources, threshold


def lead_effects(assays, leads, dbs, manifests, family, sources, threshold, hypotheses=None):
    """Extract all fixed union coordinates, with observed negatives kept intact."""
    grouped = defaultdict(list)
    for r in leads:
        grouped[r["HYPOTHESIS_ID"], r["QUERY_ID"]].append(r)

    @lru_cache(maxsize=4096)
    def at(dataset, alias, build, chrom, pos):
        key = dataset, alias
        if key not in dbs or manifests[key]["BUILD"] != build:
            return []
        return fetch_region(dbs[key], dict(CHR=chrom, START=pos, END=pos))

    effects, source_checks = [], []
    hypothesis_index = {h["HYPOTHESIS_ID"]: h for h in (hypotheses or [])}
    for (hid, qid), selected in sorted(grouped.items()):
        h = hypothesis_index.get(hid, {})
        def source_priority(q):
            key = q["DATASET"], q["ASSAY"]
            priority = 0 if key == (h.get("PRIMARY"), h.get("PRIMARY_ASSAY")) else 1 if key == (h.get("REFERENCE"), h.get("REFERENCE_ASSAY")) else 2
            return priority, q["DATASET"], q["ASSAY"], q["PROT"]
        selected.sort(key=source_priority)
        query = selected[0]
        donor, donor_assay, donor_query = None, None, None
        for q in selected:
            key = q["DATASET"], q["ASSAY"]
            rows = at(*key, q["BUILD"], q["CHR"], q["POS"])
            row, reason = source_record(rows, q["SNP"])
            if row is None and rows and not any(r["qc_reason"] == "." for r in rows):
                reason = "source_rows_failed_QC"
            m = manifests.get(key)
            if m is None:
                reason = "manifest_assay_missing"
            elif m["STATUS"] != "ready":
                reason = m["STATUS"]
            elif m["BUILD"] != q["BUILD"]:
                reason = "build_mismatch"
            agreement = compare_reported_p(q["REPORTED_P"], row, q["SELECTION_THRESHOLD"])
            source_checks.append(dict(q, FULL_SOURCE_STATUS=reason, FULL_SNP=row["snp"] if row else None,
                                      FULL_P=row["p_text"] if row else None, FULL_LOG10P=row["log10p"] if row else None,
                                      FULL_EA=row["ea"] if row else None, FULL_NEA=row["nea"] if row else None,
                                      FULL_SOURCE_ROWS_AT_COORDINATE=len(rows), FULL_SOURCE_ROW_QC_REJECTIONS=joined(r["qc_reason"] for r in rows if r["qc_reason"] != "."),
                                      SOURCE_SUMMARY_AGREEMENT=agreement,
                                      SOURCE_SELECTION_CONFIRMED=None if row is None or association_p(row) is None else int(association_p(row) <= sc.exact_p(q["SELECTION_THRESHOLD"]))))
            # Fixed source order and usable alleles, never strongest full P.
            if donor is None and row is not None and m and m["STATUS"] == "ready" and m["BUILD"] == q["BUILD"]:
                donor, donor_assay, donor_query = row, m, q
        anchor = anchor_from(donor, donor_assay) if donor else None
        for a in (r for r in assays if r["HYPOTHESIS_ID"] == hid):
            rows = at(a["DATASET"], a["ASSAY"], query["BUILD"], query["CHR"], query["POS"]) if a["MAPPED"] else []
            raw, reason = None, region_status(a, dict(BUILD=query["BUILD"]))
            aligned_beta, aligned_eaf = None, None
            if reason == "ready":
                if anchor is not None:
                    pair = sc.allele_pair_key(anchor["EA"], anchor["NEA"])
                    candidates = [r for r in rows if r["qc_reason"] == "." and sc.allele_pair_key(r["ea"], r["nea"]) == pair]
                    raw, unique_reason = strict_unique(candidates)
                    if raw is not None:
                        if (a["DATASET"], a["ASSAY"]) == (donor_assay["DATASET"], donor_assay["ASSAY"]):
                            aligned_beta, aligned_eaf, reason = raw["beta"], raw["eaf"], "anchor_source_orientation"
                        else:
                            aligned_beta, aligned_eaf, reason = si.align_record(raw, anchor, a["ANCESTRY"])
                    else:
                        reason = unique_reason if candidates else "allele_not_matched" if rows else "variant_not_observed"
                else:
                    raw, unique_reason = strict_unique([r for r in rows if r["qc_reason"] == "."])
                    reason = "source_orientation_unresolved" if raw is not None else "variant_not_observed_or_failed_QC"
            key = qid, a["PROT"], a["DATASET"], a["ASSAY"]
            selected_source = (query["BUILD"], query["CHR"], query["POS"], a["DATASET"], a["ASSAY"]) in sources
            observed_p = association_p(raw) if raw else None
            followup_pass = int(observed_p <= threshold) if anchor is not None and observed_p is not None and key in family and threshold is not None else None
            effects.append(dict(HYPOTHESIS_ID=hid, QUERY_ID=qid, REGION_ID=query["REGION_ID"], TEST_ID=uid("TEST_", key),
                                BUILD=query["BUILD"], CHR=query["CHR"], POS=query["POS"], QUERY_SNP=query["SNP"],
                                PROT=a["PROT"], DATASET=a["DATASET"], ASSAY=a["ASSAY"], ROLE=a["ROLE"], MAPPED=a["MAPPED"],
                                STATUS=reason, SOURCE_DATASET=donor_query["DATASET"] if donor_query else None,
                                SOURCE_ASSAY=donor_query["ASSAY"] if donor_query else None, SOURCE_PROT=donor_query["PROT"] if donor_query else None,
                                EA=anchor["EA"] if anchor else None, NEA=anchor["NEA"] if anchor else None,
                                OBSERVED_VARIANT=int(raw is not None), OBSERVED_HARMONIZED=int(aligned_beta is not None),
                                QUERY_VARIANT_RESOLVED=int(anchor is not None and raw is not None),
                                SOURCE_SELECTION_CONFIRMED=None if donor is None or association_p(donor) is None else int(association_p(donor) <= sc.exact_p(donor_query["SELECTION_THRESHOLD"])),
                                BETA=aligned_beta, SE=raw["se"] if raw else None, Z=aligned_beta / raw["se"] if aligned_beta is not None else None,
                                P=raw["p_text"] if raw else None, LOG10P=raw["log10p"] if raw else None, LOG10P_RAW=raw["log10p_text"] if raw else None,
                                P_STATUS=raw["p_status"] if raw else None, EAF=aligned_eaf, N=raw["n"] if raw else None,
                                RAW_SNP=raw["snp"] if raw else None, RAW_EA=raw["ea"] if raw else None, RAW_NEA=raw["nea"] if raw else None,
                                RAW_BETA=raw["beta"] if raw else None, RAW_SE=raw["se"] if raw else None, RAW_EAF=raw["eaf"] if raw else None,
                                RAW_FLAGS=raw["flags"] if raw else None,
                                SAME_SIGN_AS_SOURCE=None if aligned_beta is None or donor is None else int(aligned_beta * donor["beta"] > 0),
                                SELECTED_SOURCE=int(selected_source), FAMILY_INCLUDED=int(key in family), FAMILY_SIZE=len(family),
                                FOLLOWUP_THRESHOLD=str(threshold) if threshold else None, FOLLOWUP_PASS=followup_pass,
                                SAME_DATASET_AS_SOURCE=None if donor_assay is None else int(a["DATASET"] == donor_assay["DATASET"]),
                                INDEPENDENT_VALIDATION=0, CAPTURED_IDENTITY_CONFIRMED=0,
                                QC_REJECTIONS=joined(r["qc_reason"] for r in rows if r["qc_reason"] != ".")))
    return effects, source_checks


def write_table(path, rows, fallback):
    fields = list(rows[0]) if rows else fallback
    atomic_rows(path, rows, fields)


def run(args):
    args.out = Path(args.out)
    args.cross_dir = Path(args.cross_dir)
    args.fingerprints_dir = Path(args.fingerprints_dir) if args.fingerprints_dir else args.cross_dir / "fingerprints"
    if args.region_flank < 0 or args.cis_flank < 0 or not 0 < args.alpha < 1:
        raise ValueError("Window flanks must be nonnegative and alpha in (0,1)")
    if not args.plan_only and args.assay_manifest is None:
        raise ValueError("Full-region analysis requires --assay-manifest; --plan-only creates the frozen plan")
    args.out.mkdir(parents=True, exist_ok=True)
    # A plan in a completed output must not leave previous real-looking full
    # tables beside a new plan summary. Use a separate output or full rerun.
    summary_path = args.out / "fingerprint-check.summary.json"
    if args.plan_only and summary_path.is_file() and json.loads(summary_path.read_text()).get("mode") == "full_regions":
        raise ValueError("Use a different --out for plan-only; this directory contains completed full-region results")
    hypotheses, proteins, strict, thresholds, hashes = load_inputs(args)
    windows = sc.load_windows(args.cross_dir / "coordinates-used.tsv", args.cis_flank)
    assays, regions, leads, exclusions = build_plan(hypotheses, proteins, strict, windows, thresholds, args.region_flank)
    cis_windows = exact_cis_plan(hypotheses, windows)
    manifests, synthetic = attach_manifest(assays, args.assay_manifest)
    family, selected_sources, followup_threshold = lead_family(assays, leads, args.alpha)
    implementation = {Path(p).name: si.file_hash(p) for p in (__file__, si.__file__, sc.__file__, spec.__file__)}
    args.manifest_sha256 = si.file_hash(args.assay_manifest) if args.assay_manifest else None
    args.plan_sha256 = si.digest(dict(hypotheses=hypotheses, assays=assays, regions=regions, leads=leads,
                                    cis_windows=cis_windows, family=list(family.values()), exclusions=exclusions))
    if args.assay_manifest:
        hashes["assay_manifest"] = args.manifest_sha256
    write_table(args.out / "hypotheses.tsv", hypotheses, ["HYPOTHESIS_ID", "PROT", "ALTERNATIVE_PROT"])
    atomic_rows(args.out / "assay-plan.tsv", assays, ASSAY_FIELDS)
    atomic_rows(args.out / "regions.tsv", regions, REGION_FIELDS)
    atomic_rows(args.out / "lead-queries.tsv", leads, LEAD_FIELDS)
    write_table(args.out / "cis-windows.tsv", cis_windows, ["HYPOTHESIS_ID", "WINDOW_ID", "WINDOW_PROT", "WINDOW_ROLE", "ANNOTATION_VERIFIED"])
    write_table(args.out / "plan-exclusions.tsv", exclusions, ["HYPOTHESIS_ID", "PROT", "STATUS", "REASON"])
    write_table(args.out / "frozen-followup-family.tsv", list(family.values()), ["TEST_ID", "QUERY_ID", "PROT", "DATASET", "ASSAY", "FAMILY_SIZE"])
    requested = defaultdict(list)
    by_hypothesis = defaultdict(list)
    for r in regions:
        by_hypothesis[r["HYPOTHESIS_ID"]].append(r)
    for a in assays:
        if a["MANIFEST_STATUS"] == "ready":
            requested[a["DATASET"], a["ASSAY"]].extend((r["CHR"], r["START"], r["END"])
                for r in by_hypothesis[a["HYPOTHESIS_ID"]] if r["BUILD"] == a["BUILD"])
    dbs, receipts = {}, []
    if not args.plan_only:
        for i, (key, intervals) in enumerate(sorted(requested.items()), 1):
            if not intervals:
                continue
            print(f"fingerprint-check: {i}/{len(requested)} {key[0]}/{key[1]} complete regions", flush=True)
            db, receipt = extract_regions_cached(manifests[key], intervals, args, implementation)
            dbs[key] = db
            receipts.append(receipt)
    receipt_index = {(r["DATASET"], r["ASSAY"]): r for r in receipts}
    qc_rows, comparisons, controls, plots, cis_checks = [], [], [], [], []
    aligned_count = 0
    if not args.plan_only:
        aligned_path = args.out / "aligned-regional-variants.tsv.gz"
        aligned_fields = ["HYPOTHESIS_ID", "REGION_ID", "LEFT_PROT", "LEFT_DATASET", "LEFT_ASSAY", "RIGHT_PROT", "RIGHT_DATASET", "RIGHT_ASSAY",
                          "CHR", "POS", "EA", "NEA", "LEFT_SNP", "RIGHT_SNP", "LEFT_Z", "RIGHT_Z", "LEFT_P", "RIGHT_P", "ALIGNMENT", "IDENTIFIER_CONFLICT"]
        temporary = aligned_path.with_name("." + aligned_path.name + ".tmp")
        with gzip.open(temporary, "wt", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=aligned_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for h in hypotheses:
                hid = h["HYPOTHESIS_ID"]
                expected = [a for a in assays if a["HYPOTHESIS_ID"] == hid]
                specs, triples = comparison_specs(h, assays)
                for region in by_hypothesis[hid]:
                    @lru_cache(maxsize=3)
                    def load(key):
                        a = next(a for a in expected if assay_key(a) == key)
                        if region_status(a, region) != "ready":
                            return []
                        return fetch_region(dbs.get((a["DATASET"], a["ASSAY"])), region)
                    plotting = dict(region, assays=[])
                    for a in expected:
                        rows = load(assay_key(a))
                        qc = region_qc(hid, region, a, rows, receipt_index.get((a["DATASET"], a["ASSAY"])), thresholds[a["DATASET"]])
                        qc_rows.append(qc)
                        plotting["assays"].append(plot_profile(a, rows, qc["STATUS"]))
                    plots.append(plotting)
                    for left, right, roles in specs:
                        prefix = dict(HYPOTHESIS_ID=hid, REGION_ID=region["REGION_ID"], LEFT_PROT=left["PROT"], LEFT_DATASET=left["DATASET"], LEFT_ASSAY=left["ASSAY"],
                                      RIGHT_PROT=right["PROT"], RIGHT_DATASET=right["DATASET"], RIGHT_ASSAY=right["ASSAY"])
                        stat, aligned = compare_region(load(assay_key(left)), load(assay_key(right)), left, right)
                        statuses = ["left:" + region_status(left, region), "right:" + region_status(right, region)]
                        status = "observed_descriptive_shape" if stat["PEARSON_Z"] is not None else "insufficient_aligned_variation"
                        if any(not v.endswith(":ready") for v in statuses):
                            status = ";".join(v for v in statuses if not v.endswith(":ready"))
                        comparisons.append(dict(prefix, COMPARISON_ROLE=joined(roles), STATUS=status, **stat,
                                                ANCHOR_EXCLUDED=region["ANCHOR_EXCLUDED"], BLOCK5MB_EXCLUDED=region["BLOCK5MB_EXCLUDED"],
                                                CHROMOSOME_EXCLUDED=region["CHROMOSOME_EXCLUDED"],
                                                SAME_COHORT=int(left["DATASET"] == right["DATASET"]),
                                                SAME_PLATFORM=int(spec.PLATFORMS[left["DATASET"]] == spec.PLATFORMS[right["DATASET"]]),
                                                REGION_SELECTION="frozen_sig_union_or_gene_interval", SIGNIFICANCE_FILTER="none_within_region",
                                                COLOCALIZATION_TESTED=0, INDEPENDENT_VALIDATION=0, CAPTURED_IDENTITY_CONFIRMED=0))
                        for row in aligned:
                            writer.writerow(dict(prefix, **row))
                        aligned_count += len(aligned)
                    for primary, intended, alt in triples:
                        stat = triple_region(load(assay_key(primary)), load(assay_key(intended)), load(assay_key(alt)), primary, intended, alt)
                        state = [role + ":" + region_status(a, region) for role, a in (("primary", primary), ("intended", intended), ("alternative", alt))]
                        status = "observed_common_variant_control" if stat["DELTA_PEARSON_Z"] is not None else "insufficient_common_aligned_variation"
                        if any(not v.endswith(":ready") for v in state):
                            status = ";".join(v for v in state if not v.endswith(":ready"))
                        controls.append(dict(HYPOTHESIS_ID=hid, REGION_ID=region["REGION_ID"], PRIMARY=primary["DATASET"], PRIMARY_ASSAY=primary["ASSAY"],
                                             PROT=primary["PROT"], REFERENCE=intended["DATASET"], INTENDED_ASSAY=intended["ASSAY"], ALTERNATIVE_PROT=alt["PROT"], ALTERNATIVE_ASSAY=alt["ASSAY"],
                                             STATUS=status, **stat, ANCHOR_EXCLUDED=region["ANCHOR_EXCLUDED"], BLOCK5MB_EXCLUDED=region["BLOCK5MB_EXCLUDED"],
                                             CHROMOSOME_EXCLUDED=region["CHROMOSOME_EXCLUDED"],
                                             DESCRIPTIVE_ONLY=1, COLOCALIZATION_TESTED=0, CAPTURED_IDENTITY_CONFIRMED=0))
                    load.cache_clear()
        os.replace(temporary, aligned_path)
        effects, source_checks = lead_effects(assays, leads, dbs, manifests, family, selected_sources, followup_threshold, hypotheses)
        cis_checks = exact_cis_checks(cis_windows, assays, dbs, receipt_index, thresholds)
        write_table(args.out / "lead-effects.tsv", effects, ["HYPOTHESIS_ID", "QUERY_ID", "PROT", "DATASET", "ASSAY", "STATUS"])
        write_table(args.out / "source-summary-checks.tsv", source_checks, LEAD_FIELDS + ["FULL_SOURCE_STATUS", "SOURCE_SUMMARY_AGREEMENT"])
        write_table(args.out / "region-assay-qc.tsv", qc_rows, ["HYPOTHESIS_ID", "REGION_ID", "PROT", "DATASET", "ASSAY", "STATUS"])
        write_table(args.out / "region-comparisons.tsv", comparisons, ["HYPOTHESIS_ID", "REGION_ID", "STATUS"])
        write_table(args.out / "matched-reference-controls.tsv", controls, ["HYPOTHESIS_ID", "REGION_ID", "STATUS"])
        write_table(args.out / "cis-window-checks.tsv", cis_checks, ["HYPOTHESIS_ID", "WINDOW_ID", "PROT", "DATASET", "ASSAY", "STATUS"])
    else:
        effects, source_checks = [], []
    pairs = []
    for h in hypotheses:
        hid = h["HYPOTHESIS_ID"]
        relevant = [r for r in controls if r["HYPOTHESIS_ID"] == hid and r["REFERENCE"] == h["REFERENCE"] and r["ALTERNATIVE_ASSAY"] == h["REFERENCE_ASSAY"]]
        defined = [r for r in relevant if r["DELTA_PEARSON_Z"] is not None]
        nonanchor = [r for r in defined if not r["ANCHOR_EXCLUDED"]]
        pair_assays = [a for a in assays if a["HYPOTHESIS_ID"] == hid]
        pairs.append(dict(HYPOTHESIS_ID=hid, PROT=h["PROT"], PRIMARY=h["PRIMARY"], PRIMARY_ASSAY=h["PRIMARY_ASSAY"],
                          ALTERNATIVE_PROT=h["ALTERNATIVE_PROT"], REFERENCE=h["REFERENCE"], REFERENCE_ASSAY=h["REFERENCE_ASSAY"],
                          ORIGINAL_TIER=h["TIER"], STATUS="plan_only_full_GWAS_not_read" if args.plan_only else "descriptive_full_region_followup",
                          PLANNED_REGIONS=len(by_hypothesis[hid]), MAPPED_ASSAYS=sum(a["MAPPED"] for a in pair_assays),
                          UNMAPPED_PROTEIN_COHORTS=sum(not a["MAPPED"] for a in pair_assays),
                          FULL_ASSAYS_INSPECTED=sum((a["DATASET"], a["ASSAY"]) in dbs for a in pair_assays if a["MAPPED"]),
                          REFERENCE_CONTROL_TRIPLES=len(relevant), REFERENCE_TRIPLES_WITH_CORRELATION=len(defined),
                          REFERENCE_TRIPLES_DELTA_POSITIVE=sum(r["DELTA_PEARSON_Z"] > 0 for r in defined),
                          MEDIAN_REFERENCE_DELTA_PEARSON=float(np.median([r["DELTA_PEARSON_Z"] for r in defined])) if defined else None,
                          NONANCHOR_REFERENCE_TRIPLES=len(nonanchor), NONANCHOR_DELTA_POSITIVE=sum(r["DELTA_PEARSON_Z"] > 0 for r in nonanchor),
                          MEDIAN_NONANCHOR_DELTA_PEARSON=float(np.median([r["DELTA_PEARSON_Z"] for r in nonanchor])) if nonanchor else None,
                          SOURCE_SUMMARY_MISMATCHES=sum(r["HYPOTHESIS_ID"] == hid and r["SOURCE_SUMMARY_AGREEMENT"] in {"source_summary_P_mismatch", "source_full_P_fails_frozen_selection_threshold"} for r in source_checks),
                          SOURCE_VARIANT_IDENTITY_CONFLICTS=sum(r["HYPOTHESIS_ID"] == hid and r["FULL_SOURCE_STATUS"] == "source_variant_identifier_conflict" for r in source_checks),
                          DIRECTION_TESTED=int(any(r["HYPOTHESIS_ID"] == hid and r["N_ALIGNED"] > 0 and
                              ((r["LEFT_PROT"], r["LEFT_DATASET"], r["LEFT_ASSAY"]) == (h["PROT"], h["PRIMARY"], h["PRIMARY_ASSAY"]) or
                               (r["RIGHT_PROT"], r["RIGHT_DATASET"], r["RIGHT_ASSAY"]) == (h["PROT"], h["PRIMARY"], h["PRIMARY_ASSAY"])) for r in comparisons)),
                          COLOCALIZATION_TESTED=0, INDEPENDENT_VALIDATION=0, CAPTURED_IDENTITY_CONFIRMED=0,
                          INTERPRETATION="Descriptive, source-selected complete-region comparison; no captured-protein identity verdict"))
    write_table(args.out / "pair-summary.tsv", pairs, ["HYPOTHESIS_ID", "PROT", "ALTERNATIVE_PROT", "STATUS"])
    si.atomic_json(args.out / "input-receipts.json", receipts)
    mode = "plan_only" if args.plan_only else "full_regions"
    si.atomic_json(args.out / "plot-data.json", sc.clean_json(dict(schema="phole-fingerprint-region-plot-v1", mode=mode, synthetic_input=synthetic,
        max_points_per_assay_region=500, sampling_policy="At most 499 evenly ordered valid rows plus largest absolute Z; display only, all supplied regional rows retained in extracted gzip files", regions=plots)))
    output_names = ["hypotheses.tsv", "assay-plan.tsv", "regions.tsv", "lead-queries.tsv", "cis-windows.tsv", "plan-exclusions.tsv",
                    "frozen-followup-family.tsv", "pair-summary.tsv", "input-receipts.json", "plot-data.json"]
    if not args.plan_only:
        output_names += ["aligned-regional-variants.tsv.gz", "lead-effects.tsv", "source-summary-checks.tsv", "region-assay-qc.tsv",
                         "region-comparisons.tsv", "matched-reference-controls.tsv", "cis-window-checks.tsv"]
    outputs = [args.out / name for name in output_names]
    summary = dict(schema=SCHEMA, mode=mode, synthetic_input=synthetic, selection="explicit_pairs" if args.pairs else args.selection,
                   hypotheses=len(hypotheses), canonical_proteins=len(proteins),
                   mapped_assays=len({(a["DATASET"], a["ASSAY"]) for a in assays if a["MAPPED"]}),
                   unmapped_protein_cohorts=len({(a["PROT"], a["DATASET"]) for a in assays if not a["MAPPED"]}),
                   planned_regions=len(regions), unique_region_queries=len({(r["BUILD"], r["CHR"], r["START"], r["END"]) for r in regions}),
                   source_lead_rows=len(leads), unique_source_coordinates=len({(q["BUILD"], q["CHR"], q["POS"]) for q in leads}),
                   full_files_inspected=len(receipts), full_files_opened_this_run=sum(not r["checkpoint_reused"] for r in receipts),
                   cached_files_reused=sum(r["checkpoint_reused"] for r in receipts),
                   regional_variants_extracted=sum(r["statistics"].get("regional_rows", 0) for r in receipts),
                   aligned_regional_variant_comparisons=aligned_count, region_assay_rows=len(qc_rows), region_comparison_rows=len(comparisons),
                   matched_reference_controls=len(controls), lead_effect_rows=len(effects),
                   exact_cis_windows=len(cis_windows), exact_cis_window_assay_rows=len(cis_checks),
                   planned_unselected_coordinate_assay_tests=len(family), followup_threshold=str(followup_threshold) if followup_threshold else None,
                   source_summary_P_mismatches=sum(r["SOURCE_SUMMARY_AGREEMENT"] == "source_summary_P_mismatch" for r in source_checks),
                   source_full_selection_failures=sum(r["SOURCE_SUMMARY_AGREEMENT"] == "source_full_P_fails_frozen_selection_threshold" for r in source_checks),
                   source_variant_identity_conflicts=sum(r["FULL_SOURCE_STATUS"] == "source_variant_identifier_conflict" for r in source_checks),
                   captured_identities_confirmed=0, colocalization_tested=False,
                   complete_scope_policy="INPUT_SCOPE=full is a source declaration. Successful extraction inspects all supplied regional rows; file naming and observed variants do not prove full population variant coverage. No absence call for missing/unmapped/build-conflicting/empty/QC-invalid regions.",
                   family_definition="Frozen union source coordinate/rsID queries x every mapped intended/alternative assay; any assay selected at that coordinate excluded; missing/build-conflicting/QC-unresolved cells remain in the denominator. These are exploratory association follow-ups, not independent validation or identity tests.",
                   shape_policy="Signed Z Pearson/Spearman and sign agreement use every harmonizable allele variant without a significance filter. Matched-reference controls use an identical three-assay variant intersection. No inferential correlation P, beta-scale equivalence, LD independence, colocalization, or captured identity is asserted.",
                   holdout_policy="Region, 5 Mb, and chromosome exclusions are inherited from frozen fingerprint nominations; a region touching an exclusion is excluded in full. They are descriptive sensitivity partitions, not independent data.",
                   source_P_policy="Original text P is retained. Nonzero underflow remains precise via Decimal. Literal P=0 may be rounded, clipped or underflowed: a valid source LOG10P is used for association thresholds; zero without LOG10P has unresolved precision and is not assigned infinite strength or a threshold verdict. LOG10P is positive -log10(P); source-versus-full nonzero P agreement tolerance is 1 percent.",
                   cis_boundary_policy="Only cis-window-checks.tsv tests exact verified gene-body +/- cis_flank boundaries, including every component. Merged display-region minima must never be labelled cis minima.",
                   lead_orientation_policy="Prefer frozen nominating primary source, then nominating alternative reference, then lexical remaining source assays; choose first QC-resolved source, never strongest full P. Regional pair/triple alignment uses its own left/primary anchor and can resolve same-ancestry palindromes even if a separate lead anchor cannot.",
                   interpretation="Frozen complete-region plan only; no real complete GWAS analyzed" if args.plan_only else
                   "Descriptive follow-up of selected hypotheses. Similar fingerprints and absence of detected own-cis association do not by themselves establish that two assays captured the same protein.",
                   plan_sha256=args.plan_sha256, input_sha256=hashes, code_sha256=implementation,
                   parameters=dict(region_flank=args.region_flank, cis_flank=args.cis_flank, alpha=args.alpha, allow_scan=args.allow_scan),
                   output_sha256={p.name: si.file_hash(p) for p in sorted(outputs)}, full_gwas_identity="path/size/mtime_ns plus index fingerprint; not a full-file content hash")
    si.atomic_json(summary_path, summary)
    return summary


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cross-dir", type=Path, required=True)
    p.add_argument("--fingerprints-dir", type=Path)
    p.add_argument("--assay-manifest", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--pairs", type=Path, help="Frozen HYPOTHESIS_ID list from fingerprint assay-evidence.tsv")
    p.add_argument("--selection", choices=("core-six", "all-priority"), default="core-six",
                   help="core-six selects all current T1/T2 hypotheses (currently six; count is never hardcoded)")
    p.add_argument("--region-flank", type=int, default=250_000)
    p.add_argument("--cis-flank", type=int, default=500_000)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--allow-scan", action="store_true")
    p.add_argument("--tabix", default="tabix")
    p.add_argument("--plan-only", action="store_true")
    return p


if __name__ == "__main__":
    print(json.dumps(run(parser().parse_args()), indent=2, allow_nan=False))
