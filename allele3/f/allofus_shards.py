#!/usr/bin/env python3
"""Checkpoint one sequential ANNOVAR pass, then normalize/query independent windows.

All of Us ANNOVAR hg38 is AF-only. Windows are storage units, not independent
replicates. Missing/zero AF, unreported alleles and uncallable genotypes cannot
be converted to counts. Original gzip restart replays the prefix; no random
access claim is made for the original .idx file.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

from analysis_integrity import (InputPath, atomic_json, check_output_paths, identity,
                                transactional, validate_manifest)
from allofus_sites import (ANNOVAR_COLUMNS, ANNOVAR_POPULATIONS, BASES, Table,
                           build_name, chrom_name, frequency, integer, open_text,
                           read_targets, sort_key)

VERSION = "allofus-annovar-shards-1"
LENGTHS38 = dict(zip([str(i) for i in range(1, 23)] + ["X", "Y", "MT"],
    [248956422,242193529,198295559,190214555,181538259,170805979,159345973,
     145138636,138394717,133797422,135086622,133275309,114364328,107043718,
     101991189,90338345,83257441,80373285,58617616,64444167,46709983,
     50818468,156040895,57227415,16569]))
POPS = tuple(ANNOVAR_POPULATIONS)
AF_COLUMNS = [ANNOVAR_POPULATIONS[p]["af"] for p in POPS]
WIDE_FIELDS = ["chrom", "pos", "ref", "alt", *AF_COLUMNS, "gvs_max_af"]
LOCUS_FIELDS = ["chrom", "pos", "ref", "reported_alts", "n_reported_alts"] + [
    name + "_" + pop for pop in POPS for name in ("n_known_AF", "n_positive_ALT", "reported_AF_sum")]
TARGET_FIELDS = ["site_id", "chrom", "pos", "ref", "alt", "requested_alt", "match_status",
                 "observed_in_target", "additional_reported_alt", "source_ref", *AF_COLUMNS,
                 "gvs_max_af", "build", "evidence"]
TARGET_SITE_FIELDS = ["site_id", "chrom", "pos", "ref", "target_alts", "source_ref", "match_status",
                      "reported_alts", "additional_reported_alts", "unreported_target_alts"] + [
    name + "_" + pop for pop in POPS for name in ("n_known_AF", "n_positive_ALT")] + ["build", "evidence"]
LIMITS = [
    "AF-only; no AC, AN, callable denominator, REF count or AC-threshold event is inferred.",
    "gvs_max_af is a cross-population maximum indicator, not a ninth population.",
    "Positive AF counts reported ALT identities; AF=0 may reflect rounding and is not proof of AC=0.",
    "Missing AF, an empty source window and an unreported allele are unavailable evidence, not extinction.",
    "Marginal AF values are not renormalized or treated as a joint nucleotide probability vector.",
    "X/Y/MT retain reported frequencies; no diploidy or sample-size assumption is made.",
    "Window metrics are descriptions; gnomAD AC>=5 O/E and its denominators are not transferred.",
]


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def chromosomes(value):
    if isinstance(value, (list, tuple)):
        tokens = value
    else:
        tokens = str(value).split(",")
    selected = []
    for token in tokens:
        token = str(token).strip()
        if re.fullmatch(r"\d+-\d+", token):
            a, b = map(int, token.split("-"))
            if not 1 <= a <= b <= 22: raise ValueError("Invalid chromosome range")
            selected.extend(str(i) for i in range(a, b+1))
        else:
            selected.append(chrom_name(token))
    if not selected or len(selected) != len(set(selected)) or not set(selected) <= set(LENGTHS38):
        raise ValueError("Select distinct canonical GRCh38 chromosomes, including X/Y/MT if required")
    return sorted(selected, key=lambda c: list(LENGTHS38).index(c))


def make_plan(source, chromosomes="1-22,X,Y,MT", window_bp=10000000, source_version="unverified"):
    """Build the fixed window plan without opening or stat-ing the large source."""
    selected = globals()["chromosomes"](chromosomes)
    if isinstance(window_bp, bool) or not isinstance(window_bp, int) or window_bp < 1:
        raise ValueError("window_bp must be a positive integer")
    shards = []
    for c in selected:
        for start in range(1, LENGTHS38[c]+1, window_bp):
            end = min(start+window_bp-1, LENGTHS38[c])
            shards.append(dict(id=f"chr{c}_{start:09d}_{end:09d}", chrom=c, start=start, end=end))
    value = dict(schema=VERSION, source=os.path.abspath(os.path.expanduser(str(source))), build="GRCh38",
                 coordinate_system="1-based_inclusive", chromosomes=selected, window_bp=window_bp,
                 source_version=source_version, source_columns=ANNOVAR_COLUMNS, populations=list(POPS),
                 chromosome_lengths={c:LENGTHS38[c] for c in selected}, shards=shards)
    value["plan_id"] = json_hash(value)
    return value


def write_plan(plan, out):
    out = Path(out)
    if out.resolve() == Path(plan["source"]).resolve(): raise ValueError("Plan cannot replace source")
    if out.exists():
        old = json.loads(out.read_text())
        if old != plan: raise ValueError("Existing window plan differs; use a new plan/output directory")
        return plan
    atomic_json(out, plan)
    return plan


def read_plan(path):
    result = json.loads(Path(path).read_text())
    expected = make_plan(result["source"], result["chromosomes"], result["window_bp"], result["source_version"])
    if result != expected: raise ValueError("Window plan identity/schema or canonical GRCh38 boundaries differ")
    return result


@contextmanager
def locked_output(out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with (out.parent/("."+out.name+".lock")).open("a") as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc: raise RuntimeError("Output is in use: "+str(out)) from exc
        out.mkdir(parents=True, exist_ok=True)
        yield out


def source_rows(source):
    args = argparse.Namespace(columns=None, delimiter="tab")
    table = Table(source, args, ANNOVAR_COLUMNS)
    for number, values in table.rows():
        if table.columns != ANNOVAR_COLUMNS:
            raise ValueError("The fixed ANNOVAR preset requires the exact provided 14-column order")
        yield number, values


def raw_path(out, shard):
    return Path(out)/"raw"/(shard["id"]+".tsv.gz")


def receipt_path(out, shard):
    return Path(out)/"receipts"/(shard["id"]+".json")


def validate_receipt(out, shard, plan, source_identity):
    r = json.loads(receipt_path(out, shard).read_text())
    if (r.get("schema") != VERSION or r.get("plan_id") != plan["plan_id"] or
            r.get("shard") != shard or r.get("source_identity") != source_identity):
        raise ValueError("Prepared shard receipt has incompatible provenance: "+shard["id"])
    if r.get("status") == "empty_no_reported_rows":
        if r.get("rows") != 0 or r.get("raw_identity") is not None:
            raise ValueError("Invalid empty-window receipt")
    elif r.get("status") == "complete":
        if r.get("rows", 0) <= 0 or identity(raw_path(out, shard)) != r.get("raw_identity"):
            raise ValueError("Prepared raw shard is missing or has changed: "+shard["id"])
    else:
        raise ValueError("Prepared shard is not complete")
    return r


def prepare(plan_path, out, progress_every=1000000):
    plan = read_plan(plan_path)
    source = Path(plan["source"])
    out = check_output_paths(out, (InputPath(plan_path), InputPath(source)))
    source_id = identity(source)
    if source_id.get("status") == "missing": raise ValueError("Missing ANNOVAR source")
    if progress_every < 0: raise ValueError("progress-every must be nonnegative")
    by_window = {(s["chrom"], (s["start"]-1)//plan["window_bp"]):s for s in plan["shards"]}
    print(f"[dataset] All of Us AF-only source={source}; windows={len(plan['shards'])}; populations={','.join(POPS)}", flush=True)
    print("[build-check] All of Us preset=GRCh38; 1-based inclusive; assembly declaration, no FASTA certification", flush=True)
    with locked_output(out):
        manifest_path = out/"prepare_manifest.json"
        if manifest_path.is_file():
            old = json.loads(manifest_path.read_text())
            if old.get("plan_id") != plan["plan_id"] or old.get("source_identity") != source_id:
                raise ValueError("Prepared cache belongs to different source/plan; use a new output directory")
        else: old = {}
        (out/"raw").mkdir(exist_ok=True); (out/"receipts").mkdir(exist_ok=True)
        completed = {}
        for shard in plan["shards"]:
            if receipt_path(out, shard).is_file():
                completed[shard["id"]] = validate_receipt(out, shard, plan, source_id)
        if old.get("status") == "complete":
            if len(completed) != len(plan["shards"]): raise ValueError("Completed preparation lost a window receipt")
            print(f"[allofus-prepare] reuse complete; {len(completed)} verified windows; source gzip not reopened", flush=True)
            return old
        state = dict(schema=VERSION, status="running", plan_id=plan["plan_id"], plan=plan,
                     source_identity=source_id, limits=LIMITS,
                     original_index="ANNOVAR .idx not used; no random-access guarantee",
                     restart="Replay original gzip prefix and validate source identity; reuse completed window outputs without rewriting them")
        atomic_json(manifest_path, state)
        counters = Counter(); previous_chrom = None; previous_pos = 0; closed = set()
        current = None; handle = None; temporary = None; current_rows = 0

        def finish():
            nonlocal handle, temporary, current_rows
            if current is None: return
            sid = current["id"]
            if sid in completed:
                if completed[sid]["rows"] != current_rows:
                    raise ValueError("Replayed row count differs from completed checkpoint: "+sid)
                counters["reused_nonempty_windows"] += 1
            else:
                handle.close(); handle = None
                os.replace(temporary, raw_path(out, current)); temporary = None
                r = dict(schema=VERSION, status="complete", plan_id=plan["plan_id"], shard=current,
                         source_identity=source_id, rows=current_rows, raw_identity=identity(raw_path(out, current)))
                atomic_json(receipt_path(out, current), r); completed[sid] = r
                counters["new_nonempty_windows"] += 1
            atomic_json(out/"checkpoint.json", dict(plan_id=plan["plan_id"], completed_windows=len(completed),
                last_completed=sid, source_rows_replayed=counters["source_rows"], restart_policy=state["restart"]))

        try:
            for number, values in source_rows(source):
                counters["source_rows"] += 1
                if progress_every and counters["source_rows"] % progress_every == 0:
                    print(f"[allofus-prepare] source rows={counters['source_rows']:,}; complete windows={len(completed)}", flush=True)
                c = chrom_name(values[0])
                if c not in LENGTHS38:
                    if current is not None:
                        finish(); current = None; current_rows = 0
                    counters["noncanonical_rows"] += 1; continue
                pos = integer(values[1], f"Start on source line {number}")
                if not 1 <= pos <= LENGTHS38[c]: raise ValueError(f"Coordinate outside canonical GRCh38 {c}: {pos}")
                if c != previous_chrom:
                    if c in closed: raise ValueError("Source chromosomes must form contiguous blocks")
                    if previous_chrom is not None: closed.add(previous_chrom)
                    previous_chrom, previous_pos = c, 0
                if pos < previous_pos: raise ValueError("Source positions must be sorted within chromosome")
                previous_pos = pos
                if c not in plan["chromosomes"]:
                    if current is not None:
                        finish(); current = None; current_rows = 0
                    counters["unselected_chromosome_rows"] += 1; continue
                shard = by_window[c, (pos-1)//plan["window_bp"]]
                if current is None or current["id"] != shard["id"]:
                    finish(); current, current_rows = shard, 0
                    if shard["id"] not in completed:
                        temporary = raw_path(out, shard).with_suffix(".gz.partial")
                        handle = gzip.open(temporary, "wt", encoding="utf-8", compresslevel=1, newline="")
                        handle.write("#"+"\t".join(ANNOVAR_COLUMNS)+"\n")
                current_rows += 1; counters["selected_rows"] += 1
                if handle is not None: handle.write("\t".join(values)+"\n")
            if not counters["source_rows"]: raise ValueError("Source contains no data rows")
            finish()
            if source_id != identity(source): raise ValueError("Source changed during preparation")
            for shard in plan["shards"]:
                if shard["id"] not in completed:
                    r = dict(schema=VERSION, status="empty_no_reported_rows", plan_id=plan["plan_id"], shard=shard,
                             source_identity=source_id, rows=0, raw_identity=None)
                    atomic_json(receipt_path(out, shard), r); completed[shard["id"]] = r
            state.update(status="complete", counts=dict(counters), n_windows=len(completed),
                         n_nonempty_windows=sum(r["rows"] > 0 for r in completed.values()))
            atomic_json(manifest_path, state)
            print(f"[allofus-prepare] complete rows={counters['source_rows']:,}; nonempty windows={state['n_nonempty_windows']}; empty windows mean no reported source rows", flush=True)
            return state
        except BaseException as exc:
            if handle is not None: handle.close()
            state.update(status="failed", error=str(exc), counts=dict(counters), completed_windows=len(completed))
            atomic_json(manifest_path, state)
            raise


def prepared_info(prepared_dir, shard_id=None):
    directory = Path(prepared_dir)
    m = json.loads((directory/"prepare_manifest.json").read_text())
    if m.get("status") != "complete" or m.get("schema") != VERSION:
        raise ValueError("Preparation must complete before independent normalization/query")
    plan = m["plan"]
    if plan != make_plan(plan["source"], plan["chromosomes"], plan["window_bp"], plan["source_version"]):
        raise ValueError("Preparation plan identity has changed")
    if m.get("plan_id") != plan["plan_id"]: raise ValueError("Preparation identity mismatch")
    if shard_id is None: return m, None, None
    selected = [s for s in plan["shards"] if s["id"] == shard_id]
    if len(selected) != 1: raise ValueError("Unknown shard ID")
    s = selected[0]
    return m, s, validate_receipt(directory, s, plan, m["source_identity"])


def value_text(value):
    return "NA" if value is None else str(value)


class BGZFWriter:
    def __init__(self, path, fields):
        import pysam
        self.path, self.fields = Path(path), fields
        self.handle = pysam.BGZFile(str(path), "wb")
        self.handle.write(("#"+"\t".join(fields)+"\n").encode())
        self.rows = 0
    def write(self, row):
        self.handle.write(("\t".join(value_text(row.get(k)) for k in self.fields)+"\n").encode())
        self.rows += 1
    def close(self):
        self.handle.close()
        if self.rows:
            import pysam
            pysam.tabix_index(str(self.path), seq_col=0, start_col=1, end_col=1, meta_char="#", zerobased=False, force=True)


def parse_snv(values, line_number):
    c = chrom_name(values[0]); pos = integer(values[1], "SNV Start")
    ref, alt = values[3].upper(), values[4].upper()
    if ref not in BASES or alt not in BASES: return None
    if ref == alt: raise ValueError(f"REF equals ALT on line {line_number}")
    if integer(values[2], "SNV End") != pos: raise ValueError("SNV Start and inclusive End differ")
    metrics = {col:frequency(values[ANNOVAR_COLUMNS.index(col)], col) for col in [*AF_COLUMNS, "gvs_max_af"]}
    return dict(chrom=c, pos=pos, ref=ref, alt=alt, **metrics)


def normalize_shard(prepared_dir, shard_id, out):
    m, shard, receipt = prepared_info(prepared_dir, shard_id)
    path = Path(out)/"manifest.json"
    if path.is_file():
        old = json.loads(path.read_text())
        if old.get("status") == "complete_descriptive_unverified":
            old = validate_manifest(out, require_evidence=True)
            if old.get("plan_id") != m["plan_id"] or old.get("shard") != shard or old.get("raw_identity") != receipt["raw_identity"]:
                raise ValueError("Normalized output belongs to a different prepared shard")
            print(f"[allofus-normalize] reuse {shard_id}", flush=True)
            return old
    return _normalize_shard(prepared_dir, shard_id, out)


@transactional("manifest.json", input_dirs=("prepared_dir",), input_prefixes=())
def _normalize_shard(prepared_dir, shard_id, out):
    m, shard, receipt = prepared_info(prepared_dir, shard_id)
    out = Path(out)
    print(f"[dataset] All of Us shard={shard_id}; source_rows={receipt['rows']}; AF-only", flush=True)
    print("[build-check] GRCh38 exact prepared window; 1-based SNVs; X/Y/MT remain frequency summaries", flush=True)
    aw = BGZFWriter(out/"wide_alleles.tsv.gz", WIDE_FIELDS)
    sw = BGZFWriter(out/"wide_loci.tsv.gz", LOCUS_FIELDS)
    counts = Counter(); per_pop = {p:Counter() for p in POPS}; locus = {}; current = None; previous = 0
    def emit():
        if not locus: return
        rows = sorted(locus.values(), key=lambda r:r["alt"])
        counts["snv_loci"] += 1
        summary = {k:rows[0][k] for k in ("chrom", "pos", "ref")}
        summary.update(reported_alts=",".join(r["alt"] for r in rows), n_reported_alts=len(rows))
        for pop, col in zip(POPS, AF_COLUMNS):
            known = [r[col] for r in rows if r[col] is not None]
            positive = sum(x > 0 for x in known)
            summary.update({"n_known_AF_"+pop:len(known), "n_positive_ALT_"+pop:positive,
                            "reported_AF_sum_"+pop:sum(known) if known else None})
            per_pop[pop]["loci"] += 1
            per_pop[pop]["loci_all_reported_AF_known"] += len(known) == len(rows)
            per_pop[pop]["loci_ge2_positive_ALT"] += positive >= 2
            per_pop[pop]["loci_ge3_positive_ALT"] += positive >= 3
            per_pop[pop]["loci_AF_sum_gt1"] += bool(known) and sum(known) > 1 + len(rows)*1e-6
            per_pop[pop]["positive_ALT_rows"] += positive
        for row in rows: aw.write(row)
        sw.write(summary)
    try:
        if receipt["rows"]:
            for number, values in source_rows(raw_path(prepared_dir, shard)):
                counts["source_rows"] += 1
                c = chrom_name(values[0]); pos = integer(values[1], "Start")
                if c != shard["chrom"] or not shard["start"] <= pos <= shard["end"] or pos < previous:
                    raise ValueError("Raw row violates its shard coordinate/order contract")
                previous = pos
                row = parse_snv(values, number)
                if row is None: counts["non_snv_rows"] += 1; continue
                key = c, pos
                if current != key:
                    emit(); locus = {}; current = key
                elif locus and next(iter(locus.values()))["ref"] != row["ref"]:
                    raise ValueError("Conflicting canonical SNV REF at the same coordinate")
                if row["alt"] in locus:
                    if locus[row["alt"]] != row: raise ValueError("Conflicting duplicate SNV frequency metrics")
                    counts["identical_duplicate_rows"] += 1
                else:
                    locus[row["alt"]] = row; counts["unique_snv_alleles"] += 1
            emit()
        if counts["source_rows"] != receipt["rows"]: raise ValueError("Raw shard row count differs from receipt")
    finally:
        aw.close(); sw.close()
    return dict(status="complete_descriptive_unverified", schema=VERSION, plan_id=m["plan_id"], shard=shard,
                source_identity=m["source_identity"], raw_identity=receipt["raw_identity"], counts=dict(counts),
                data_status="reported_snvs" if counts["unique_snv_alleles"] else "no_reported_snv_rows",
                populations={p:dict(v) for p,v in per_pop.items()}, build="GRCh38", evidence="frequency_only",
                cache_schema="one wide row per distinct SNV allele; one wide row per locus; never population-expanded genome-wide output",
                index_status="BGZF_tabix_1based" if counts["unique_snv_alleles"] else "no_index_for_empty_table",
                limits=LIMITS)


def read_normalized(normalized_root, shard, plan_id, source_identity):
    directory = Path(normalized_root)/shard["id"]
    m = validate_manifest(directory, require_evidence=True)
    if (m.get("schema") != VERSION or m.get("plan_id") != plan_id or m.get("shard") != shard
            or m.get("source_identity") != source_identity):
        raise ValueError("Normalized shard provenance differs: "+shard["id"])
    return directory, m


@transactional("manifest.json", input_files=("plan_path",), input_dirs=("normalized_root",), input_prefixes=())
def summarize(plan_path, normalized_root, out):
    plan = read_plan(plan_path); source_identity = None; totals = Counter(); pop_total = {p:Counter() for p in POPS}
    rows = []; receipts = []
    for shard in plan["shards"]:
        directory = Path(normalized_root)/shard["id"]
        m = validate_manifest(directory, require_evidence=True)
        if source_identity is None: source_identity = m.get("source_identity")
        directory, m = read_normalized(normalized_root, shard, plan["plan_id"], source_identity)
        totals.update(m["counts"])
        receipts.append(identity(directory/"manifest.json", full=True))
        for pop in POPS:
            pop_total[pop].update(m["populations"][pop])
            rows.append(dict(shard_id=shard["id"], chrom=shard["chrom"], start=shard["start"], end=shard["end"],
                             population=pop, data_status=m["data_status"], **m["populations"][pop]))
    fields = ["shard_id", "chrom", "start", "end", "population", "data_status", "loci", "loci_all_reported_AF_known",
              "loci_ge2_positive_ALT", "loci_ge3_positive_ALT", "loci_AF_sum_gt1", "positive_ALT_rows"]
    write_tsv(Path(out)/"population_window_summary.tsv", fields, rows)
    write_tsv(Path(out)/"population_summary.tsv", ["population", *fields[6:]],
              [dict(population=p, **v) for p,v in pop_total.items()])
    return dict(status="complete_descriptive_unverified", schema=VERSION, plan_id=plan["plan_id"],
                build="GRCh38", evidence="frequency_only", source_identity=source_identity,
                counts=dict(totals), n_windows=len(plan["shards"]), input_manifests=receipts,
                populations={p:dict(v) for p,v in pop_total.items()}, limits=LIMITS)


def write_tsv(path, fields, rows):
    op = gzip.open(path, "wt", newline="", compresslevel=3) if str(path).endswith(".gz") else Path(path).open("w", newline="")
    with op as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows: writer.writerow({k:value_text(row.get(k)) for k in fields})


@transactional("manifest.json", input_dirs=("prepared_dir", "normalized_root"), input_files=("targets_path",), input_prefixes=())
def query(prepared_dir, normalized_root, targets_path, out):
    m, _, _ = prepared_info(prepared_dir)
    plan = m["plan"]; targets = read_targets(targets_path, "GRCh38", "GRCh38")
    windows = {(s["chrom"], (s["start"]-1)//plan["window_bp"]):s for s in plan["shards"]}
    grouped = {}; unavailable = []; counts = Counter(); output = []; site_output = []; used = []
    for key,target in targets.items():
        shard = windows.get((key[0], (key[1]-1)//plan["window_bp"]))
        if shard is None or key[1] > shard["end"]:
            unavailable.append((target, "outside_prepared_windows"))
        else: grouped.setdefault(shard["id"], (shard, []))[1].append(target)

    def emit(target, records, state="queried"):
        requested = set(target["alts"].split(",")); source_refs = {r["ref"] for r in records}
        if len(source_refs) > 1: raise ValueError("Normalized cache has conflicting SNV references")
        source_ref = next(iter(source_refs), None)
        compatible = source_ref == target["ref"]
        if source_ref is not None and not compatible: state = "reference_mismatch"
        elif not records and state == "queried": state = "not_in_summary"
        reported = {r["alt"]:r for r in records} if compatible else {}
        all_alts = sorted(requested | set(reported))
        for alt in all_alts:
            r = reported.get(alt)
            known_metric = r is not None and any(r[col] is not None for col in AF_COLUMNS)
            status = ("matched_frequency_only" if known_metric else "matched_metrics_missing") if r else state if state != "queried" else "not_in_summary"
            output.append(dict(site_id=target["site_id"], chrom=target["chrom"], pos=target["pos"], ref=target["ref"],
                alt=alt, requested_alt=int(alt in requested), observed_in_target=int(alt in requested),
                additional_reported_alt=int(alt not in requested and r is not None), match_status=status,
                source_ref=source_ref, **({k:r[k] for k in [*AF_COLUMNS,"gvs_max_af"]} if r else {}),
                build="GRCh38", evidence="frequency_only" if known_metric else "metrics_missing" if r else "unavailable"))
            counts[status] += 1
        summary = dict(site_id=target["site_id"], chrom=target["chrom"], pos=target["pos"], ref=target["ref"],
            target_alts=target["alts"], source_ref=source_ref,
            match_status="matched_coordinate" if compatible else state,
            reported_alts=",".join(sorted(reported)) if compatible else None,
            additional_reported_alts=",".join(sorted(set(reported)-requested)) if compatible else None,
            unreported_target_alts=",".join(sorted(requested-set(reported))), build="GRCh38",
            evidence=("frequency_only" if any(r[col] is not None for r in reported.values() for col in AF_COLUMNS)
                      else "metrics_missing") if compatible else "unavailable")
        for pop,col in zip(POPS, AF_COLUMNS):
            summary["n_known_AF_"+pop] = sum(r[col] is not None for r in reported.values()) if compatible else None
            summary["n_positive_ALT_"+pop] = sum(r[col] is not None and r[col] > 0 for r in reported.values()) if compatible else None
        site_output.append(summary)

    import pysam
    for sid,(shard,items) in grouped.items():
        directory, nm = read_normalized(normalized_root, shard, plan["plan_id"], m["source_identity"])
        used.append(identity(directory/"manifest.json", full=True))
        if not nm["counts"].get("unique_snv_alleles", 0):
            for t in items: emit(t, [], "no_reported_snv_rows_in_window")
            continue
        path = directory/"wide_alleles.tsv.gz"
        with pysam.TabixFile(str(path)) as tabix:
            for target in sorted(items, key=lambda t:t["pos"]):
                records = []
                for line in tabix.fetch(target["chrom"], target["pos"]-1, target["pos"]):
                    values = line.split("\t")
                    if len(values) != len(WIDE_FIELDS): raise ValueError("Invalid normalized cache schema")
                    row = dict(zip(WIDE_FIELDS,values)); row["pos"] = integer(row["pos"], "cache position")
                    if row["pos"] != target["pos"]: continue
                    for col in [*AF_COLUMNS,"gvs_max_af"]: row[col] = frequency(row[col],col)
                    records.append(row)
                emit(target, records)
    for target,state in unavailable: emit(target, [], state)
    output.sort(key=lambda r:(*sort_key((r["chrom"],r["pos"])),r["alt"]))
    site_output.sort(key=lambda r:sort_key((r["chrom"],r["pos"])))
    write_tsv(Path(out)/"target_alleles.tsv.gz",TARGET_FIELDS,output)
    write_tsv(Path(out)/"target_sites.tsv",TARGET_SITE_FIELDS,site_output)
    return dict(status="complete_descriptive_unverified", schema=VERSION, plan_id=plan["plan_id"],
                build="GRCh38", evidence="frequency_only", target_count=len(targets), queried_windows=len(grouped),
                source_identity=m["source_identity"], source_manifest=identity(Path(prepared_dir)/"prepare_manifest.json",full=True),
                normalized_manifests=used, counts=dict(counts), targets=identity(targets_path,full=True), limits=LIMITS,
                target_flags="requested_alt/observed_in_target indicate the supplied candidate ALT list, not individual genotypes; additional_reported_alt means source-catalog ALT outside that list, not a new mutation",
                scan_policy="Only candidate windows and indexed exact positions; original billion-row gzip not opened")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("plan", help="Create canonical GRCh38 windows without opening the source")
    q.add_argument("--source", required=True); q.add_argument("--out", required=True)
    q.add_argument("--chromosomes", default="1-22,X,Y,MT"); q.add_argument("--window-bp", type=int, default=10000000)
    q.add_argument("--source-version", default="unverified")
    q = sub.add_parser("prepare", help="One sequential gzip pass; persistent per-window checkpoints")
    q.add_argument("--plan", required=True); q.add_argument("--out", required=True)
    q.add_argument("--progress-every", type=int, default=1000000)
    q = sub.add_parser("normalize", help="Normalize a single completed raw window to indexed wide tables")
    q.add_argument("--prepared", required=True); q.add_argument("--shard", required=True); q.add_argument("--out", required=True)
    q = sub.add_parser("summarize", help="Require all normalized planned windows and aggregate small receipts")
    q.add_argument("--plan", required=True); q.add_argument("--normalized-root", required=True); q.add_argument("--out", required=True)
    q = sub.add_parser("query", help="Read only candidate positions from completed indexed wide windows")
    q.add_argument("--prepared", required=True); q.add_argument("--normalized-root", required=True)
    q.add_argument("--targets", required=True); q.add_argument("--out", required=True)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command == "plan":
        result = write_plan(make_plan(args.source,args.chromosomes,args.window_bp,args.source_version),args.out)
        print(json.dumps(dict(status="planned",plan_id=result["plan_id"],windows=len(result["shards"]))))
    elif args.command == "prepare": result = prepare(args.plan,args.out,args.progress_every)
    elif args.command == "normalize": result = normalize_shard(args.prepared,args.shard,args.out)
    elif args.command == "summarize": result = summarize(args.plan,args.normalized_root,args.out)
    else: result = query(args.prepared,args.normalized_root,args.targets,args.out)
    if args.command != "plan": print(json.dumps(dict(status=result["status"],counts=result.get("counts",{})),sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
