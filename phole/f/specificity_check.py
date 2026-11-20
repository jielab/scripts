#!/usr/bin/env python3
"""Check externally selected cis variants in full GWAS without inferring assay identity.

The external reference fixes the variant and its true EA/NEA columns. Every
manifest assay of the same canonical protein is inspected, including observed
nonsignificant effects. A missing or unorientable variant is never a zero.
The frozen Bonferroni family counts planned, unselected variant--assay tests,
including those that later prove missing or fail QC. It is an exploratory
follow-up family, not a test that an antibody captured another protein.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

import signed_input as si

SCHEMA = "phole-specificity-full-gwas-v1"
STANDARD = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
KNOWN_TARGETS = {"single_gene", "complex_complete"}


def finite(value):
    return value is not None and isinstance(value, (float, int)) and math.isfinite(value)


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def exact_p(text):
    try:
        p = Decimal(str(text))
    except InvalidOperation:
        raise ValueError("invalid_P") from None
    if not p.is_finite() or not 0 <= p <= 1:
        raise ValueError("invalid_P")
    return p


def dataset(value):
    value = str(value).lower().strip()
    value = "ukb" if value == "prot" else value
    if value not in si.DATASETS:
        raise ValueError(f"Unknown dataset {value!r}")
    return value


def boolean(value, field):
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ValueError(f"{field} must be a definite true/false or 1/0")


def read_queries(path, query_build=None):
    output, seen = [], set()
    for index, raw in enumerate(si.read_tsv(path), 1):
        r = {k: str(v).strip() for k, v in raw.items() if k is not None}
        r["QUERY_ELIGIBLE"] = boolean(r.get("QUERY_ELIGIBLE", "1"), "QUERY_ELIGIBLE")
        if not r["QUERY_ELIGIBLE"]:
            r.setdefault("QUERY_ID", f"excluded_row_{index}")
            output.append(r)
            continue
        for col in ("QUERY_ID", "PROT", "REFERENCE_ASSAY", "CHR", "POS", "SNP"):
            if not r.get(col) or r[col] == ".":
                raise ValueError(f"Query requires {col}: {r.get('QUERY_ID', '?')}")
        if r["QUERY_ID"] in seen:
            raise ValueError(f"Duplicate QUERY_ID: {r['QUERY_ID']}")
        seen.add(r["QUERY_ID"])
        r["REFERENCE_DATASET"] = dataset(r.get("REFERENCE_DATASET", r.get("REFERENCE", "")))
        r["PRIMARY_DATASET"] = dataset(r.get("PRIMARY_DATASET", r.get("PRIMARY", "ukb")))
        r["REFERENCE_ASSAY"] = si.safe_assay(r["REFERENCE_ASSAY"])
        if r.get("PRIMARY_ASSAY") not in (None, "", "."):
            r["PRIMARY_ASSAY"] = si.safe_assay(r["PRIMARY_ASSAY"])
        r["CHR"] = si.normal_chr(r["CHR"])
        try:
            r["POS"] = int(r["POS"])
        except ValueError:
            raise ValueError(f"Query {r['QUERY_ID']}: POS must be a positive integer") from None
        if r["CHR"] not in STANDARD or r["POS"] <= 0:
            raise ValueError(f"Query {r['QUERY_ID']}: invalid coordinate")
        r["BUILD"] = si.normal_build(r.get("BUILD", query_build))
        output.append(r)
    return output


def load_windows(path, flank):
    """Read effective annotation, not SNP query coordinates; expand once."""
    grouped = defaultdict(list)
    if path is None:
        return {}
    for r in si.read_tsv(path):
        p = r.get("prot", r.get("PROT", ""))
        if not p:
            raise ValueError("Cis annotation requires prot or PROT")
        grouped[p].append(r)
    result = {}
    for p, rows in grouped.items():
        valid, reasons, build = [], [], None
        for r in rows:
            status = r.get("annotation_status", r.get("STATUS", ""))
            if status not in KNOWN_TARGETS:
                reasons.append("target_identity_or_components_unresolved")
            if "identity_verified" in r:
                try:
                    if not boolean(r["identity_verified"], "identity_verified"):
                        reasons.append("target_identity_not_verified")
                except ValueError:
                    reasons.append("target_identity_verification_ambiguous")
            try:
                b = si.normal_build(r.get("build", r.get("BUILD", "")))
                c = si.normal_chr(r.get("CHR", ""))
                lo, hi = int(r["START"]), int(r["END"])
                if c not in STANDARD or lo < 1 or hi < lo:
                    raise ValueError
                if build is not None and b != build:
                    reasons.append("mixed_annotation_builds")
                build = b
                valid.append((c, max(1, lo - flank), hi + flank))
            except (ValueError, KeyError):
                reasons.append("target_coordinate_unresolved")
        # Merge overlapping intervals to avoid duplicate tabix output/counting.
        intervals = []
        for c, lo, hi in sorted(set(valid), key=lambda v: (si.chrom_order(v[0]), v[1], v[2])):
            if intervals and intervals[-1][0] == c and lo <= intervals[-1][2] + 1:
                intervals[-1] = (c, intervals[-1][1], max(hi, intervals[-1][2]))
            else:
                intervals.append((c, lo, hi))
        result[p] = dict(intervals=intervals, build=build, verified=bool(intervals) and not reasons,
                         reason=";".join(sorted(set(reasons))) or "verified_effective_target_intervals")
    return result


def strict_unique(records):
    row, status = si.unique_record(records)
    if row is None:
        return row, status
    base = exact_p(row["p_text"])
    # Shared input utility compares numeric fields with a small absolute
    # tolerance. Tiny P values require a relative-only duplicate check here.
    for other in records[1:]:
        p = exact_p(other["p_text"])
        if abs(p - base) > max(p, base) * Decimal("1e-7"):
            return None, "duplicate_precise_P_conflict"
    return row, status


def reference_record(records, requested_snp):
    chosen = [r for r in records if r["snp"] == requested_snp or si.snp_tokens(r["snp"]) & si.snp_tokens(requested_snp)]
    if chosen:
        return strict_unique(chosen)
    if len(records) == 1:
        if si.snp_tokens(requested_snp) and si.snp_tokens(records[0]["snp"]):
            return None, "reference_variant_identifier_conflict"
        return records[0], "unique_coordinate_ID_differs"
    return strict_unique(records)


def oriented_record(records, anchor, ancestry):
    aligned, rejected = [], Counter()
    for row in records:
        beta, eaf, reason = si.align_record(row, anchor, ancestry)
        if beta is None:
            rejected[reason] += 1
        else:
            aligned.append(dict(row, beta=beta, eaf=eaf, ea=anchor["EA"], nea=anchor["NEA"], alignment=reason))
    if not aligned:
        return None, next(iter(rejected), "variant_absent_or_failed_row_QC"), dict(rejected)
    row, status = strict_unique(aligned)
    return row, row["alignment"] if row is not None else status, dict(rejected)


def raw_at_anchor(records, anchor):
    """Association P does not require a known strand orientation.

    A palindromic SNP can have an identifiable allele pair and a readable P
    while its signed effect cannot be aligned across ancestries. Preserve both
    facts, without assigning an oriented beta to that cell.
    """
    pair = allele_pair_key(anchor["EA"], anchor["NEA"])
    compatible = [r for r in records if allele_pair_key(r["ea"], r["nea"]) == pair]
    return strict_unique(compatible)


def allele_pair_key(ea, nea):
    pair = tuple(sorted((ea, nea)))
    if len(ea) == len(nea) == 1:
        pair = min(pair, tuple(sorted((ea.translate(si.COMPLEMENT), nea.translate(si.COMPLEMENT)))))
    return pair


def empty_window(reason="not_requested"):
    return dict(requested=False, input_complete=False, variant_rows=0, valid_P_rows=0, valid_effect_rows=0,
                distinct_positions=0, min_P=None, min_P_SNP=None, min_P_CHR=None, min_P_POS=None,
                intervals=[], intervals_with_rows=0, reason=reason, has_gws=None,
                interpretation="No claim about absence of a cis association")


def extract_with_windows(assay, positions, window, allow_scan, tabix):
    """Inspect full declared GWAS windows and exact positions in a single pass."""
    path = Path(assay["PATH"])
    header = si.read_header(path)
    records, stats, bad_positions = defaultdict(list), Counter(), defaultdict(Counter)
    intervals = window["intervals"]
    hits = Counter()
    distinct_positions = set()
    best = None
    window_result = empty_window()
    window_result.update(requested=True, intervals=intervals, reason=window["reason"])
    column_chr, column_pos, column_p = header.index("CHR"), header.index("POS"), header.index("P")

    def consume(line):
        nonlocal best
        if not line.strip() or line.startswith("#"):
            return
        fields = line.split()
        if len(fields) >= 2 and fields[0].upper() == header[0] and fields[1].upper() == header[1]:
            return
        stats["rows_read"] += 1
        try:
            if len(fields) != len(header):
                raise ValueError("malformed_columns")
            c, pos = si.normal_chr(fields[column_chr]), int(fields[column_pos])
            if c not in STANDARD or pos < 1:
                raise ValueError("invalid_coordinate")
        except (ValueError, IndexError):
            stats["unlocatable_rows"] += 1
            return
        indexes = [i for i, (ch, lo, hi) in enumerate(intervals) if ch == c and lo <= pos <= hi]
        requested = (c, pos) in positions
        if not requested and not indexes:
            return
        if indexes:
            stats["window_variant_rows"] += 1
            distinct_positions.add((c, pos))
            for i in indexes:
                hits[i] += 1
            try:
                p = exact_p(fields[column_p])
                stats["window_valid_P_rows"] += 1
                if best is None or p < best[0]:
                    snp = fields[header.index("SNP")] if "SNP" in header else "."
                    best = p, fields[column_p], snp, c, pos
            except ValueError:
                stats["window_invalid_P_rows"] += 1
        try:
            row = si.parse_record(fields, header, assay.get("N"))
            exact_p(row["p_text"])
        except ValueError as e:
            stats["rejected_" + str(e)] += 1
            if requested:
                bad_positions[c, pos][str(e)] += 1
            if indexes:
                stats["window_row_QC_rejected"] += 1
            return
        if indexes:
            stats["window_valid_effect_rows"] += 1
        if requested:
            records[c, pos].append(row)

    index_paths = [Path(assay["INDEX"])] if assay.get("INDEX") else [Path(str(path) + ".tbi"), Path(str(path) + ".csi")]
    index = next((p for p in index_paths if p.is_file()), None)
    executable = shutil.which(tabix)
    if index is not None and executable is not None:
        indexed_path = str(path) if index in (Path(str(path) + ".tbi"), Path(str(path) + ".csi")) else str(path) + "##idx##" + str(index)
        listed = subprocess.run([executable, "-l", indexed_path], text=True, capture_output=True, check=True).stdout.splitlines()
        contigs = defaultdict(list)
        for c in listed:
            contigs[si.normal_chr(c)].append(c)
        if any(len(v) > 1 for v in contigs.values()):
            raise ValueError(f"{path}: duplicated normalized tabix contigs")
        regions = list(intervals) + [(c, p, p) for c, p in positions if not any(ch == c and lo <= p <= hi for ch, lo, hi in intervals)]
        with tempfile.TemporaryDirectory(prefix="phole-specificity-") as temporary:
            reg_path, err_path = Path(temporary) / "regions.tsv", Path(temporary) / "stderr"
            with reg_path.open("w") as handle:
                for c, lo, hi in sorted(regions, key=lambda v: (si.chrom_order(v[0]), v[1], v[2])):
                    if c in contigs:
                        handle.write(f"{contigs[c][0]}\t{lo}\t{hi}\n")
                    else:
                        stats["requested_intervals_absent_contig"] += 1
            if reg_path.stat().st_size:
                with err_path.open("w") as error:
                    process = subprocess.Popen([executable, "-R", str(reg_path), indexed_path], stdout=subprocess.PIPE, stderr=error, text=True)
                    try:
                        for line in process.stdout:
                            consume(line)
                    finally:
                        process.stdout.close()
                        code = process.wait()
                if code:
                    raise RuntimeError(f"tabix failed for {path}: {err_path.read_text()[-2000:]}")
        method = "tabix_exact_positions_and_complete_cis_intervals"
    else:
        if not allow_scan:
            raise ValueError(f"{path}: full cis window scan requires a usable tabix index or explicit --allow-scan")
        with si.open_text(path) as handle:
            for line in handle:
                consume(line)
        method = "explicit_full_file_scan_with_cis_intervals"
    window_result.update(variant_rows=stats["window_variant_rows"], valid_P_rows=stats["window_valid_P_rows"],
                         valid_effect_rows=stats["window_valid_effect_rows"], distinct_positions=len(distinct_positions),
                         intervals_with_rows=len(hits), min_P=best[1] if best else None, min_P_SNP=best[2] if best else None,
                         min_P_CHR=best[3] if best else None, min_P_POS=best[4] if best else None,
                         has_gws=bool(best[0] <= Decimal("5e-8")) if best else None)
    window_result["input_complete"] = bool(window["verified"] and intervals and len(hits) == len(intervals)
                                           and not stats["unlocatable_rows"] and not stats["window_row_QC_rejected"])
    window_result["interpretation"] = ("Significant variant observed in queried cis region" if window_result["has_gws"] else
        "No P <= 5e-8 among variants in the fully inspected declared GWAS cis intervals; not proof of no biological effect"
        if window_result["input_complete"] else "Cis absence unresolved: incomplete annotation/variant coverage/QC")
    if stats["unlocatable_rows"]:
        window_result["reason"] += ";unlocatable_rows_in_read_stream"
    if stats["window_row_QC_rejected"]:
        window_result["reason"] += ";cis_rows_failed_QC"
    if len(hits) != len(intervals):
        window_result["reason"] += ";one_or_more_intervals_without_variants"
    return records, {"method": method, "header": header, "index": str(index) if index else None, **dict(stats)}, window_result, bad_positions


def extract_cached(assay, positions, window, args, implementation):
    fingerprint = si.summary_identity(assay)
    specification = dict(schema=SCHEMA, assay=assay, positions=sorted(positions), window=window,
                         input=fingerprint, implementation=implementation, scan_window=args.scan_cis_window)
    key = si.digest(specification)
    checkpoint = (getattr(args, "state_dir", None) or args.out / "state") / assay["DATASET"] / (si.digest(assay["ASSAY"])[:24] + ".json")
    if checkpoint.is_file():
        saved = json.loads(checkpoint.read_text())
        if saved.get("key") == key:
            records = defaultdict(list)
            for item in saved["records"]:
                for row in item["rows"]:
                    for field in ("n", "eaf"):
                        row[field] = math.nan if row[field] is None else row[field]
                records[item["chr"], item["pos"]] = item["rows"]
            bad = {(r["chr"], r["pos"]): r["reasons"] for r in saved["bad_positions"]}
            receipt = saved["receipt"] | {"checkpoint_reused": True}
            return records, saved["window"], bad, receipt
    if args.scan_cis_window and window and window["build"] == assay["BUILD"] and window["intervals"]:
        records, stats, w, bad = extract_with_windows(assay, positions, window, args.allow_scan, args.tabix)
    else:
        records, stats = si.extract_positions(assay, positions, args.allow_scan, args.tabix)
        reason = "not_requested" if not args.scan_cis_window else "annotation_missing_or_build_incompatible"
        w, bad = empty_window(reason), {}
    if si.summary_identity(assay) != fingerprint:
        raise ValueError(f"Input changed during extraction: {assay['PATH']}")
    receipt = dict(DATASET=assay["DATASET"], ASSAY=assay["ASSAY"], PROT=assay["PROT"], BUILD=assay["BUILD"],
                   fingerprint=fingerprint, statistics=stats, checkpoint_key=key, checkpoint_reused=False)
    encoded = [dict(chr=c, pos=p, rows=rr) for (c, p), rr in sorted(records.items())]
    bad_encoded = [dict(chr=c, pos=p, reasons=dict(v)) for (c, p), v in sorted(bad.items())]
    si.atomic_json(checkpoint, clean_json(dict(key=key, records=encoded, bad_positions=bad_encoded, window=w, receipt=receipt)))
    return records, w, bad, receipt


def run(args):
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = si.validate_manifest(args.assay_manifest)
    all_queries = read_queries(args.queries, getattr(args, "query_build", None))
    primary = getattr(args, "primary", "ukb")
    if primary not in (*si.DATASETS, "all"):
        raise ValueError("primary must be ukb/ckb/decode/scot/all")
    excluded_queries = [q for q in all_queries if not q["QUERY_ELIGIBLE"]]
    queries = [q for q in all_queries if q["QUERY_ELIGIBLE"] and (primary == "all" or q["PRIMARY_DATASET"] == primary)]
    synthetic = {r["SYNTHETIC"] for r in manifest}
    if len(synthetic) > 1:
        raise ValueError("Synthetic and real GWAS inputs cannot be mixed")
    if args.scan_cis_window and args.cis_annotation is None:
        raise ValueError("--scan-cis-window requires --cis-annotation with effective gene intervals")
    if args.cis_flank < 0 or not 0 < args.alpha < 1:
        raise ValueError("cis-flank must be nonnegative and alpha between 0 and 1")
    assays = {(r["DATASET"], r["ASSAY"]): r for r in manifest}
    by_protein = defaultdict(list)
    for row in manifest:
        by_protein[row["PROT"]].append(row)
    requested = defaultdict(set)
    query_issues = {}
    for q in queries:
        donor_key = q["REFERENCE_DATASET"], q["REFERENCE_ASSAY"]
        donor = assays.get(donor_key)
        issues = []
        if donor is None:
            issues.append("reference_assay_not_in_manifest")
        elif donor["PROT"] != q["PROT"]:
            raise ValueError(f"Query {q['QUERY_ID']}: reference assay belongs to a different PROT")
        elif donor["BUILD"] != q["BUILD"]:
            issues.append("reference_query_build_mismatch")
        if q.get("PRIMARY_ASSAY") not in (None, "", "."):
            primary_assay = assays.get((q["PRIMARY_DATASET"], q["PRIMARY_ASSAY"]))
            if primary_assay is None:
                issues.append("primary_assay_not_in_manifest")
            elif primary_assay["PROT"] != q["PROT"]:
                raise ValueError(f"Query {q['QUERY_ID']}: primary assay belongs to a different PROT")
        query_issues[q["QUERY_ID"]] = issues
        for assay in by_protein[q["PROT"]]:
            if assay["BUILD"] == q["BUILD"]:
                requested[assay["DATASET"], assay["ASSAY"]].add((q["CHR"], q["POS"]))
    windows = load_windows(args.cis_annotation, args.cis_flank)
    implementation = {Path(__file__).name: si.file_hash(__file__), Path(si.__file__).name: si.file_hash(si.__file__)}
    extracted, window_results, rejected_positions, receipts = {}, {}, {}, []

    def extract_one(key):
        assay = assays[key]
        records, w, bad, receipt = extract_cached(assay, requested[key], windows.get(assay["PROT"]), args, implementation)
        extracted[key], window_results[key], rejected_positions[key] = records, w, bad
        receipts.append(receipt)

    donor_keys = {(q["REFERENCE_DATASET"], q["REFERENCE_ASSAY"]) for q in queries} & set(requested)
    for key in sorted(donor_keys):
        extract_one(key)
    anchors, variants, donor_records = {}, {}, {}
    for q in queries:
        qid, dk = q["QUERY_ID"], (q["REFERENCE_DATASET"], q["REFERENCE_ASSAY"])
        rows = extracted.get(dk, {}).get((q["CHR"], q["POS"]), [])
        row, reason = reference_record(rows, q["SNP"])
        if any(i.startswith("reference_") for i in query_issues[qid]):
            row, reason = None, ";".join(query_issues[qid])
        if row is None:
            anchors[qid] = None
            variants[qid] = (q["PROT"], q["BUILD"], q["CHR"], q["POS"], "unresolved", dk, q["SNP"])
        else:
            donor = assays[dk]
            anchors[qid] = dict(CHR=q["CHR"], POS=q["POS"], EA=row["ea"], NEA=row["nea"], EAF=row["eaf"],
                                ANCESTRY=donor["ANCESTRY"], BUILD=q["BUILD"])
            variants[qid] = (q["PROT"], q["BUILD"], q["CHR"], q["POS"], allele_pair_key(row["ea"], row["nea"]))
            try:
                reference_threshold = exact_p(q.get("REFERENCE_THRESHOLD", "5e-8"))
                if not Decimal(0) < reference_threshold <= Decimal("5e-8"):
                    raise ValueError("reference threshold must be positive and no larger than 5e-8")
                if exact_p(row["p_text"]) > reference_threshold:
                    query_issues[qid].append("reference_full_P_fails_original_selection_threshold")
            except ValueError:
                raise ValueError(f"Query {qid}: invalid REFERENCE_THRESHOLD") from None
        donor_records[qid] = row, reason
    # Everything in this set was selected to nominate a query, and cannot also
    # be counted as external validation in another duplicate query row.
    selected = {(variants[q["QUERY_ID"]], q["REFERENCE_DATASET"], q["REFERENCE_ASSAY"]) for q in queries}
    family = set()
    for q in queries:
        for assay in by_protein[q["PROT"]]:
            key = variants[q["QUERY_ID"]], assay["DATASET"], assay["ASSAY"]
            if key not in selected:
                family.add(key)
    threshold = Decimal(str(args.alpha)) / len(family) if family else None
    frozen = [dict(TEST_ID=si.digest(key)[:24], PROT=key[0][0], BUILD=key[0][1], CHR=key[0][2], POS=key[0][3],
                   VARIANT_KEY=json.dumps(key[0]), DATASET=key[1], ASSAY=key[2], FAMILY_SIZE=len(family),
                   FOLLOWUP_THRESHOLD=str(threshold)) for key in sorted(family, key=lambda v: json.dumps(v))]
    si.atomic_tsv(args.out / "frozen-followup-family.tsv", frozen,
                  list(frozen[0]) if frozen else ["TEST_ID", "PROT", "DATASET", "ASSAY", "FAMILY_SIZE"])
    # The family depends only on frozen queries, reference variant identities
    # and the planned assay manifest. Recipient P values never enter it.
    for key in sorted(set(requested) - donor_keys):
        extract_one(key)
    effects = []
    query_summaries = []
    for q in queries:
        qid, variant, anchor = q["QUERY_ID"], variants[q["QUERY_ID"]], anchors[q["QUERY_ID"]]
        donor, donor_reason = donor_records[qid]
        donor_selected = bool(donor and exact_p(donor["p_text"]) <= exact_p(q.get("REFERENCE_THRESHOLD", "5e-8")))
        target_window = windows.get(q["PROT"])
        cis_verified = (bool(target_window["verified"] and target_window["build"] == q["BUILD"] and
                             any(c == q["CHR"] and lo <= q["POS"] <= hi for c, lo, hi in target_window["intervals"]))
                        if target_window else None)
        primary_rows = []
        for assay in sorted(by_protein[q["PROT"]], key=lambda a: (a["DATASET"], a["ASSAY"])):
            key = assay["DATASET"], assay["ASSAY"]
            comparison = variant, *key
            reference_selected = comparison in selected
            w = window_results.get(key, empty_window("assay_build_incompatible"))
            record, raw_record, orientation, rejects = None, None, "unresolved_reference_variant", {}
            if assay["BUILD"] != q["BUILD"]:
                orientation = "recipient_query_build_mismatch"
            elif anchor is not None:
                raw_records = extracted.get(key, {}).get((q["CHR"], q["POS"]), [])
                raw_record, _ = raw_at_anchor(raw_records, anchor)
                if key == (q["REFERENCE_DATASET"], q["REFERENCE_ASSAY"]):
                    record, raw_record, orientation = donor, donor, "selected_reference_EA"
                else:
                    record, orientation, rejects = oriented_record(raw_records, anchor, assay["ANCESTRY"])
            association_record = record if record is not None else raw_record
            row_qc = rejected_positions.get(key, {}).get((q["CHR"], q["POS"]), {})
            if record is None and row_qc:
                rejects = {**rejects, **row_qc}
            same = bool(record["beta"] * donor["beta"] > 0) if record is not None and donor is not None else None
            passed = bool(exact_p(association_record["p_text"]) <= threshold) if association_record is not None and threshold is not None and not reference_selected else None
            if reference_selected:
                state = ("selected_reference_not_validation" if donor_selected else "donor_full_association_not_confirmed") if record is not None else "reference_unresolved"
            elif association_record is None:
                state = "untestable"
            elif not donor_selected:
                state = "donor_full_association_not_confirmed"
            elif passed and record is None:
                state = "association_direction_unresolved"
            elif passed:
                state = "association_same_direction" if same else "association_opposite_or_zero_direction"
            else:
                state = "not_significant_at_frozen_followup_threshold"
            r = dict(QUERY_ID=qid, TEST_ID=si.digest(comparison)[:24], PROT=q["PROT"],
                     PRIMARY_DATASET=q["PRIMARY_DATASET"], PRIMARY_ASSAY=q.get("PRIMARY_ASSAY", "."),
                     REFERENCE_DATASET=q["REFERENCE_DATASET"], REFERENCE_ASSAY=q["REFERENCE_ASSAY"],
                     REFERENCE_QUERY_SNP=q["SNP"], REFERENCE_FULL_SNP=donor["snp"] if donor else ".",
                     REFERENCE_STATUS=donor_reason, REFERENCE_P=donor["p_text"] if donor else ".",
                     REFERENCE_REPORTED_P=q.get("REFERENCE_P", "."), REFERENCE_SELECTION_THRESHOLD=q.get("REFERENCE_THRESHOLD", "5e-8"),
                     REFERENCE_STILL_SELECTED="." if donor is None else int(donor_selected),
                     REFERENCE_CIS_ANNOTATION_VERIFIED="." if cis_verified is None else int(cis_verified),
                     REFERENCE_BETA=donor["beta"] if donor else ".", REFERENCE_SE=donor["se"] if donor else ".",
                     DATASET=assay["DATASET"], ASSAY=assay["ASSAY"], BUILD=q["BUILD"], CHR=q["CHR"], POS=q["POS"],
                     SNP=association_record["snp"] if association_record else ".", EA=anchor["EA"] if anchor else ".", NEA=anchor["NEA"] if anchor else ".",
                     OBSERVED_VARIANT=int(association_record is not None),
                     OBSERVED_HARMONIZED=int(record is not None), ALIGNMENT=orientation,
                     BETA=record["beta"] if record else ".", SE=association_record["se"] if association_record else ".",
                     Z=record["beta"] / record["se"] if record else ".", P=association_record["p_text"] if association_record else ".",
                     EAF=record["eaf"] if record and finite(record["eaf"]) else ".", N=record["n"] if record and finite(record["n"]) else ".",
                     RAW_EA=raw_record["ea"] if raw_record else ".", RAW_NEA=raw_record["nea"] if raw_record else ".",
                     RAW_BETA=raw_record["beta"] if raw_record else ".", RAW_SE=raw_record["se"] if raw_record else ".",
                     RAW_P=raw_record["p_text"] if raw_record else ".",
                     RAW_EAF=raw_record["eaf"] if raw_record and finite(raw_record["eaf"]) else ".",
                     RAW_N=raw_record["n"] if raw_record and finite(raw_record["n"]) else ".",
                     EFFECT_SCALE=assay.get("EFFECT_SCALE", "unknown"), ANCESTRY=assay["ANCESTRY"],
                     SAME_DIRECTION="." if same is None else int(same), FAMILY_INCLUDED=int(comparison in family),
                     FAMILY_SIZE=len(family), FOLLOWUP_THRESHOLD=str(threshold) if threshold is not None else ".",
                     FOLLOWUP_PASS="." if passed is None else int(passed), STATUS=state,
                     CIS_REPLICATION_SUPPORTED=int(bool(donor_selected and cis_verified and passed and same and
                         not reference_selected and assay["DATASET"] != q["REFERENCE_DATASET"])),
                     REFERENCE_SELECTED_IN_QUERY_LIST=int(reference_selected),
                     SAME_DATASET_AS_REFERENCE=int(assay["DATASET"] == q["REFERENCE_DATASET"]),
                     CIS_WINDOW_REQUESTED=int(w["requested"]), CIS_WINDOW_COMPLETE=int(w["input_complete"]),
                     CIS_WINDOW_VARIANT_ROWS=w["variant_rows"], CIS_WINDOW_VALID_P_ROWS=w["valid_P_rows"],
                     CIS_WINDOW_MIN_P=w["min_P"] or ".", CIS_WINDOW_MIN_SNP=w["min_P_SNP"] or ".",
                     CIS_WINDOW_HAS_GWS="." if w["has_gws"] is None else int(w["has_gws"]),
                     CIS_WINDOW_INTERPRETATION=w["interpretation"], CIS_WINDOW_REASON=w["reason"],
                     QUERY_ISSUES=";".join(query_issues[qid]) or ".", ROW_REJECTIONS=json.dumps(rejects, sort_keys=True),
                     HETEROGENEITY_TEST="not_performed_effect_scales_and_sampling_covariance_unverified",
                     SYNTHETIC=next(iter(synthetic)), WRONG_PROTEIN_VERDICT="not_established")
            effects.append(r)
            if assay["DATASET"] == q["PRIMARY_DATASET"] and q.get("PRIMARY_ASSAY") in (None, "", ".", assay["ASSAY"]):
                primary_rows.append(r)
        query_summaries.append(dict(QUERY_ID=qid, PROT=q["PROT"], PRIMARY=q["PRIMARY_DATASET"],
                                   PRIMARY_ASSAY=q.get("PRIMARY_ASSAY", "."), REFERENCE=q["REFERENCE_DATASET"],
                                   REFERENCE_ASSAY=q["REFERENCE_ASSAY"], CHR=q["CHR"], POS=q["POS"], SNP=q["SNP"],
                                   REFERENCE_RESOLVED=int(anchor is not None), PRIMARY_ASSAYS_INSPECTED=len(primary_rows),
                                   REFERENCE_FULL_SELECTION_CONFIRMED=int(donor_selected),
                                   REFERENCE_CIS_ANNOTATION_VERIFIED="." if cis_verified is None else int(cis_verified),
                                   PRIMARY_ASSAYS_OBSERVED=sum(r["OBSERVED_VARIANT"] for r in primary_rows),
                                   PRIMARY_ASSAYS_HARMONIZED=sum(r["OBSERVED_HARMONIZED"] for r in primary_rows),
                                   PRIMARY_FOLLOWUP_ASSOCIATIONS=sum(r["FOLLOWUP_PASS"] == 1 for r in primary_rows),
                                   PRIMARY_OPPOSITE_ASSOCIATIONS=sum(r["FOLLOWUP_PASS"] == 1 and r["SAME_DIRECTION"] == 0 for r in primary_rows),
                                   PRIMARY_FULL_CIS_WINDOWS=sum(r["CIS_WINDOW_COMPLETE"] for r in primary_rows),
                                   PRIMARY_CIS_SIGNIFICANT_WINDOWS=sum(r["CIS_WINDOW_HAS_GWS"] == 1 for r in primary_rows),
                                   QUERY_ISSUES=";".join(query_issues[qid]) or ".", WRONG_PROTEIN_VERDICT="not_established"))
    window_table = []
    for key, w in sorted(window_results.items()):
        assay = assays[key]
        window_table.append(dict(DATASET=key[0], ASSAY=key[1], PROT=assay["PROT"], BUILD=assay["BUILD"],
                                 **{k: json.dumps(v) if isinstance(v, list) else "." if v is None else v for k,v in w.items()}))
    si.atomic_tsv(args.out / "targeted-cis-effects.tsv", effects, list(effects[0]) if effects else ["QUERY_ID", "PROT", "DATASET", "ASSAY", "STATUS"])
    si.atomic_tsv(args.out / "specificity-query-results.tsv", query_summaries, list(query_summaries[0]) if query_summaries else ["QUERY_ID", "PROT"])
    si.atomic_tsv(args.out / "cis-window-checks.tsv", window_table, list(window_table[0]) if window_table else ["DATASET", "ASSAY", "PROT"])
    si.atomic_tsv(args.out / "specificity-excluded-queries.tsv",
                  [dict(QUERY_ID=q.get("QUERY_ID", "."), PROT=q.get("PROT", "."), REASON="QUERY_ELIGIBLE=0") for q in excluded_queries],
                  ["QUERY_ID", "PROT", "REASON"])
    si.atomic_json(args.out / "input-receipts.json", clean_json(receipts))
    summary = dict(schema=SCHEMA, synthetic=next(iter(synthetic)) == "TRUE", primary=primary,
                   input_query_rows=len(all_queries), ineligible_query_rows=len(excluded_queries), queries=len(queries), canonical_proteins=len({q["PROT"] for q in queries}),
                   reference_queries_resolved=sum(a is not None for a in anchors.values()), planned_unselected_variant_assay_tests=len(family),
                   reference_full_selection_confirmed=sum(r["REFERENCE_FULL_SELECTION_CONFIRMED"] for r in query_summaries),
                   threshold=str(threshold) if threshold is not None else None, output_rows=len(effects),
                   unique_observed_unselected_tests=len({r["TEST_ID"] for r in effects if r["FAMILY_INCLUDED"] and r["OBSERVED_VARIANT"]}),
                   unique_harmonized_unselected_tests=len({r["TEST_ID"] for r in effects if r["FAMILY_INCLUDED"] and r["OBSERVED_HARMONIZED"]}),
                   unique_significant_unselected_tests=len({r["TEST_ID"] for r in effects if r["FOLLOWUP_PASS"] == 1}),
                   unique_cis_replication_supported_tests=len({r["TEST_ID"] for r in effects if r["CIS_REPLICATION_SUPPORTED"]}),
                   fully_inspected_cis_assays=sum(w["input_complete"] for w in window_results.values()),
                   wrong_proteins_confirmed=0,
                   interpretation="Targeted full-GWAS follow-up only. Neither a significance difference nor an absent cis effect establishes a wrong captured protein.",
                   family_definition="Deduplicated canonical protein/build/coordinate/true allele pair x recipient assay; reference-selected comparisons excluded; missing, wrong-build and unorientable planned tests retained. Unresolved donor tests remain conservative planned tests. Frozen before recipient P inspection.",
                   cis_window_scope="Only declared full-summary variants within effective annotated windows. Variant coverage, ancestry and assay reliability remain distinct from absence of biological effect.",
                   effect_scale_policy="No cross-platform beta difference, heterogeneity or equivalence test without verified comparable phenotype scales and covariance.",
                   hashes={"queries_sha256":si.file_hash(args.queries), "assay_manifest_sha256":si.file_hash(args.assay_manifest),
                           "cis_annotation_sha256":si.file_hash(args.cis_annotation) if args.cis_annotation else None, **implementation},
                   full_gwas_identity="path/size/mtime_ns and index fingerprints; not a full GWAS content hash",
                   parameters={"primary":primary,"alpha":args.alpha,"cis_flank":args.cis_flank,"scan_cis_window":args.scan_cis_window,"allow_scan":args.allow_scan})
    si.atomic_json(args.out / "specificity-check.summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assay-manifest", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, help="Reusable specificity checkpoints")
    parser.add_argument("--primary", choices=(*si.DATASETS, "all"), default="ukb")
    parser.add_argument("--query-build", choices=("37", "38"))
    parser.add_argument("--allow-scan", action="store_true")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--cis-annotation", type=Path)
    parser.add_argument("--scan-cis-window", action="store_true")
    parser.add_argument("--cis-flank", type=int, default=500_000)
    parser.add_argument("--alpha", type=float, default=0.05)
    print(json.dumps(run(parser.parse_args()), indent=2))


if __name__ == "__main__":
    main()
