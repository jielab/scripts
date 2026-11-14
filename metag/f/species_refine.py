#!/usr/bin/env python3
"""Second-stage species rescue by competitive short-read alignment.

Why this stage exists
---------------------
Kraken is excellent for broad read classification, but a read ending at genus /
family / order (or taxid 0 under a nonzero confidence threshold) should not be
*forced* to species merely because many genomes exist in the reference library.
This program uses species already supported by the sample as a candidate set,
selects one deterministic representative assembly per species, and competitively
aligns unresolved/higher-rank paired reads with minimap2 ``-x sr``.

A rescue is accepted only when the primary alignment has sufficient MAPQ,
aligned fraction and identity.  For a Kraken higher-rank call, the candidate
species must also descend from that original taxon.  Conflicting mates are not
rescued.  Original species calls are never overwritten.

Rescues are candidate-panel-dependent evidence, not validated species truth.
MAPQ measures competition only among the selected references; omitted close
relatives can give a misleadingly unique alignment. Candidate limits and skipped
genomes are therefore reported explicitly. Novel/unrepresented reads can still
map to related references, so an increased assigned fraction is not accuracy.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import json
import math
import mmap
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import threading
import time

from fastq_stream import fragments, read_id
from result_output import managed_outputs, prepare as prepare_result, finish as finish_result, replace_bool
from common import atomic_json, check_space, file_state, load_config, require_compute, sha256
from bin_layout import fastq_name
from species_split import BufferedWriters, _cigar_metrics
from taxonomy import Taxonomy, safe_name
from assignment_io import open_assignment_text

_U64 = struct.Struct("<Q")
_TAXON_RE = re.compile(r"^T(\d+)\|")


def _executable(value):
    path = Path(value).expanduser() if os.sep in value else None
    resolved = str(path.resolve(strict=True)) if path is not None else shutil.which(value)
    if not resolved or not Path(resolved).is_file() or not os.access(resolved, os.X_OK):
        raise RuntimeError(f"Cannot find executable: {value}")
    return str(Path(resolved).resolve())


def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", newline="")


def _open_fasta(path):
    path = Path(path)
    return gzip.open(path, "rb") if path.name.endswith(".gz") else path.open("rb")


def _load_representatives(path, reference_raw):
    result = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species_taxid", "name", "group", "path", "bases"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise RuntimeError(f"Representative table missing {required}: {path}")
        for row in reader:
            tid = int(row["species_taxid"])
            bases = int(row["bases"])
            if tid <= 0 or bases <= 0:
                raise RuntimeError("Representative taxids and base counts must be positive")
            raw = Path(row["path"])
            raw = raw if raw.is_absolute() else Path(reference_raw) / raw
            raw = raw.resolve(strict=True)
            if tid in result:
                raise RuntimeError(f"Duplicate species representative taxid {tid}")
            if not raw.is_file():
                raise RuntimeError(f"Representative is not a regular file: {raw}")
            result[tid] = {**row, "species_taxid": tid, "bases": bases, "resolved_path": str(raw)}
    return result


def _load_candidate_counts(path):
    candidates = Counter()
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species_taxid", "status", "fragments"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise RuntimeError(f"counts.tsv missing {required}: {path}")
        for row in reader:
            if row["status"] != "assigned_species":
                continue
            tid = int(row["species_taxid"] or 0)
            if tid:
                candidates[tid] += int(row["fragments"])
    return candidates



def _extra_taxids_from_sylph(path, min_ani=95.0):
    if not path:
        return []
    result = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species_taxid", "adjusted_ani"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise RuntimeError(f"Sylph candidate table missing {required}: {path}")
        for row in reader:
            try:
                ani = float(row["adjusted_ani"])
                if math.isfinite(ani) and min_ani <= ani <= 100:
                    result.append(int(row["species_taxid"]))
            except (TypeError, ValueError):
                continue
    return sorted(set(result))

def choose_candidates(counts, representatives, *, minimum, maximum, max_reference_gib,
                      max_single_reference_gib, extra=()):
    ranked = [(count, tid) for tid, count in counts.items() if count >= minimum and tid in representatives]
    for tid in extra:
        if tid in representatives:
            ranked.append((max(minimum, int(counts.get(tid, 0))), tid))
    ranked.sort(key=lambda x: (-x[0], x[1]))
    selected = []
    skipped = [{"species_taxid": tid, "reason": "representative_missing", "bases": 0,
                "kraken_fragments": int(counts.get(tid, 0))}
               for tid in sorted({tid for tid, n in counts.items() if n >= minimum} | set(extra))
               if tid not in representatives]
    budget = int(max_reference_gib * 2**30)
    single = int(max_single_reference_gib * 2**30)
    used = 0
    seen = set()
    for count, tid in ranked:
        if tid in seen:
            continue
        seen.add(tid)
        row = representatives[tid]
        bases = int(row["bases"])
        if bases > single:
            skipped.append({"species_taxid": tid, "reason": "single_reference_too_large", "bases": bases, "kraken_fragments": count})
            continue
        if len(selected) >= maximum:
            skipped.append({"species_taxid": tid, "reason": "max_species", "bases": bases, "kraken_fragments": count})
            continue
        if used + bases > budget:
            skipped.append({"species_taxid": tid, "reason": "reference_budget", "bases": bases, "kraken_fragments": count})
            continue
        selected.append({**row, "kraken_fragments": int(counts.get(tid, 0))})
        used += bases
    if not selected:
        raise RuntimeError("No candidate species fit the requested evidence/reference limits")
    return selected, skipped, used


def filter_human_candidates(counts, extra, representatives, *, include_human):
    """Apply the human exclusion to every candidate source, including explicit IDs."""
    if include_human:
        return counts, list(extra)
    def human(tid):
        return tid == 9606 or representatives.get(tid, {}).get("group") == "human"
    return Counter({tid: n for tid, n in counts.items() if not human(tid)}), [tid for tid in extra if not human(tid)]


def build_competitive_reference(selected, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    mapping = []
    with destination.open("wb", buffering=8 << 20) as out:
        for row in selected:
            tid = int(row["species_taxid"])
            records = bases = 0
            with _open_fasta(row["resolved_path"]) as inp:
                for line in inp:
                    if line.startswith(b">"):
                        records += 1
                        out.write(f">T{tid}|R{records}\n".encode())
                    else:
                        seq = line.rstrip(b"\r\n")
                        if seq:
                            bases += len(seq)
                            out.write(seq + b"\n")
            if not records or not bases:
                raise RuntimeError(f"Empty candidate reference: {row['resolved_path']}")
            mapping.append({"species_taxid": tid, "name": row["name"], "group": row["group"],
                            "source_path": row["resolved_path"], "records": records, "bases_written": bases,
                            "kraken_fragments": row["kraken_fragments"]})
    return mapping


def _stop(process):
    if process is None or process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)


def _sam_positive(fields, *, min_mapq, min_fraction, min_identity):
    flag = int(fields[1])
    if flag & 0x4:
        return False, 0
    mapq = int(fields[4])
    aligned, consumed = _cigar_metrics(fields[5])
    sequence = fields[9]
    read_length = len(sequence) if sequence != "*" else consumed
    if read_length <= 0 or aligned <= 0:
        return False, 0
    nm = None
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            nm = int(field[5:]); break
    if nm is None:
        raise RuntimeError("minimap2 primary SAM row lacks NM:i tag")
    identity = max(0.0, 1.0 - nm / aligned)
    fraction = aligned / read_length
    match = _TAXON_RE.match(fields[2])
    if match is None:
        raise RuntimeError(f"Competitive reference contig lacks T<taxid>| prefix: {fields[2]}")
    tid = int(match.group(1))
    return (mapq >= min_mapq and fraction >= min_fraction and identity >= min_identity), tid


def align_for_rescue(r1, r2, assignments, reference, destination, *, expected,
                     minimap2, threads, rescue_statuses, min_mapq, min_fraction, min_identity,
                     index_batch_gib=32):
    """Stream only ambiguous fragments to minimap2; store one candidate taxid per original ordinal."""
    destination = Path(destination)
    calls_path, eligible_path, seen_path = (destination / name for name in
        ("rescue-taxids.bin", "eligible.bin", "seen.bin"))
    for p in (calls_path, eligible_path, seen_path):
        if p.exists() or p.is_symlink():
            raise RuntimeError(f"Refusing to overwrite refinement state: {p}")
    calls_size = expected * 8
    with calls_path.open("w+b") as cf, eligible_path.open("w+b") as ef, seen_path.open("w+b") as sf:
        cf.truncate(calls_size); ef.truncate(expected); sf.truncate(expected)
        calls = mmap.mmap(cf.fileno(), calls_size)
        eligible = mmap.mmap(ef.fileno(), expected)
        seen = mmap.mmap(sf.fileno(), expected)
        errors, result = [], {}
        process = producer = None
        if index_batch_gib <= 0:
            raise ValueError("index_batch_gib must be positive")
        command = [minimap2, "-ax", "sr", "--secondary=no", "-I", f"{index_batch_gib}G",
                   "-t", str(threads), str(reference), "-"]
        log_path = destination / "species-refine.minimap2.log"
        started = time.monotonic()

        def producer_fn():
            count = eligible_count = 0
            stream = fragments(r1, r2)
            try:
                with open_assignment_text(assignments) as ah:
                    reader = csv.DictReader(ah, delimiter="\t")
                    required = {"ordinal", "r1_read_id", "r2_read_id", "status", "assigned_taxid"}
                    if not reader.fieldnames or not required <= set(reader.fieldnames):
                        raise RuntimeError(f"Assignments missing {required}")
                    status_col = "kraken_status" if "kraken_status" in reader.fieldnames else "status"
                    for count, pair in enumerate(stream, 1):
                        row = next(reader, None)
                        if row is None or int(row["ordinal"]) != count:
                            raise RuntimeError(f"FASTQ/assignment ordinal mismatch at {count}")
                        a, b = pair
                        if row["r1_read_id"] != read_id(a[0]).decode("utf-8", "backslashreplace"):
                            raise RuntimeError(f"R1 ID differs at ordinal {count}")
                        if b is not None and row["r2_read_id"] != read_id(b[0]).decode("utf-8", "backslashreplace"):
                            raise RuntimeError(f"R2 ID differs at ordinal {count}")
                        if row["status"] == "host_rescued" or row[status_col] not in rescue_statuses:
                            continue
                        eligible[count - 1] = 1
                        eligible_count += 1
                        for mate, record in enumerate((a, b), 1):
                            if record is None:
                                continue
                            header = f"@{count}/{mate}\n".encode()
                            process.stdin.write(header + record[1].rstrip(b"\r\n") + b"\n+\n" +
                                                record[3].rstrip(b"\r\n") + b"\n")
                        if eligible_count % 500_000 == 0:
                            print(f"species rescue streamed {eligible_count:,} eligible / {count:,} total fragments", flush=True)
                    if count != expected:
                        raise RuntimeError(f"FASTQ fragments differ from expected: {count} != {expected}")
                    if next(reader, None) is not None:
                        raise RuntimeError("Assignments contain more rows than FASTQ")
                result.update(total_fragments=count, eligible_fragments=eligible_count)
                process.stdin.close()
            except BaseException as exc:
                errors.append(exc); _stop(process)
                try: process.stdin.close()
                except (OSError, ValueError): pass
            finally:
                stream.close()

        try:
            with log_path.open("w", encoding="utf-8") as log:
                log.write("COMMAND " + json.dumps(command) + "\n"); log.flush()
                process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                           stderr=log, start_new_session=True, bufsize=1 << 20)
                producer = threading.Thread(target=producer_fn, daemon=True, name="metag-species-refine-producer")
                producer.start()
                primary = positive = conflicts = 0
                for raw in process.stdout:
                    if raw.startswith(b"@"):
                        continue
                    fields_b = raw.rstrip(b"\r\n").split(b"\t")
                    if len(fields_b) < 11:
                        raise RuntimeError(f"Malformed minimap2 SAM row: {raw[:160]!r}")
                    fields = [x.decode("utf-8", "strict") for x in fields_b]
                    if int(fields[1]) & (0x100 | 0x800):
                        continue
                    qname = fields[0]
                    ordinal_text = qname.split("/", 1)[0]
                    if not ordinal_text.isdigit():
                        raise RuntimeError(f"Unexpected minimap2 qname: {qname}")
                    ordinal = int(ordinal_text)
                    if not 1 <= ordinal <= expected or not eligible[ordinal - 1]:
                        raise RuntimeError(f"SAM ordinal not eligible/out of range: {ordinal}")
                    flag = int(fields[1])
                    bit = 1 if flag & 0x40 else (2 if flag & 0x80 else 0)
                    if not bit:
                        # minimap2 may omit paired flags if it did not recognise the pair;
                        # the /1,/2 suffix still gives an unambiguous mate identity.
                        bit = 1 if qname.endswith("/1") else (2 if qname.endswith("/2") else 0)
                    if not bit:
                        raise RuntimeError(f"Cannot determine mate for SAM row {qname}")
                    idx = ordinal - 1
                    if seen[idx] & bit:
                        raise RuntimeError(f"Duplicate primary SAM mate at ordinal {ordinal}")
                    seen[idx] |= bit; primary += 1
                    ok, tid = _sam_positive(fields, min_mapq=min_mapq,
                                            min_fraction=min_fraction, min_identity=min_identity)
                    if not ok:
                        continue
                    positive += 1
                    offset = idx * 8
                    current = _U64.unpack_from(calls, offset)[0]
                    if seen[idx] & 0x4:
                        continue
                    if current == 0 or current == tid:
                        _U64.pack_into(calls, offset, tid)
                    else:
                        _U64.pack_into(calls, offset, 0)
                        seen[idx] |= 0x4
                        conflicts += 1
                code = process.wait(); producer.join(timeout=30)
                if producer.is_alive():
                    raise RuntimeError("Species-refine FASTQ producer did not stop")
                if errors:
                    raise RuntimeError(f"Species-refine producer failed: {errors[0]}") from errors[0]
                if code:
                    raise RuntimeError(f"minimap2 exited {code}; see {log_path}")
                missing = 0
                for i in range(expected):
                    if eligible[i] and (seen[i] & 0x3) != 0x3:
                        missing += 1
                if missing:
                    raise RuntimeError(f"minimap2 primary SAM coverage incomplete for {missing:,} eligible fragments")
                calls.flush(); cf.flush(); os.fsync(cf.fileno())
                result.update(primary_alignments=primary, positive_alignments=positive,
                              conflicting_fragments=conflicts, seconds=time.monotonic() - started,
                              calls=str(calls_path), calls_sha256=sha256(calls_path))
                atomic_json(destination / "alignment-summary.json", result)
        finally:
            _stop(process)
            if producer is not None: producer.join(timeout=10)
            if process is not None:
                if process.stdout: process.stdout.close()
                try: process.stdin.close()
                except (OSError, ValueError, AttributeError): pass
            calls.close(); eligible.close(); seen.close()
    eligible_path.unlink(missing_ok=True); seen_path.unlink(missing_ok=True)
    return result


def _call_at(mapping, index):
    return _U64.unpack_from(mapping, index * 8)[0]


def apply_rescue(assignments, calls_path, taxonomy, expected, destination):
    """Apply candidate calls with lineage guard and write a refined assignment table."""
    counts = Counter(); metadata = {}; rejected_lineage = accepted = positive = 0
    if Path(calls_path).stat().st_size != expected * 8:
        raise RuntimeError("Refinement call count differs from the expected fragment count")
    with Path(calls_path).open("rb") as raw:
        calls = mmap.mmap(raw.fileno(), expected * 8, access=mmap.ACCESS_READ)
        try:
            with open_assignment_text(assignments) as inp, gzip.open(destination, "wt", compresslevel=3,
                                                            encoding="utf-8", newline="") as out:
                reader = csv.DictReader(inp, delimiter="\t")
                if not reader.fieldnames:
                    raise RuntimeError("Assignments have no header")
                status_col = "kraken_status" if "kraken_status" in reader.fieldnames else "status"
                tid_col = "kraken_taxid" if "kraken_taxid" in reader.fieldnames else "assigned_taxid"
                fields = ["ordinal", "r1_read_id", "r2_read_id", "assigned_taxid", "assigned_rank",
                          "species_taxid", "name", "group", "status", "source", "original_taxid",
                          "original_status", "alignment_candidate_taxid", "lineage_consistent"]
                writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                n = 0
                for n, row in enumerate(reader, 1):
                    if int(row["ordinal"]) != n or n > expected:
                        raise RuntimeError(f"Assignment ordinal mismatch at {n}")
                    original_tid = int(row[tid_col]); original_status = row[status_col]
                    candidate = _call_at(calls, n - 1)
                    if candidate:
                        candidate_info = taxonomy.describe(candidate)
                        if candidate_info["status"] != "assigned_species" or candidate_info["species_taxid"] != candidate:
                            raise RuntimeError(f"Alignment candidate is not an exact species taxid: {candidate}")
                    positive += int(bool(candidate))
                    consistent = True
                    final_tid, source = original_tid, "kraken"
                    if row.get("status") == "host_rescued":
                        final_tid, source = int(row["assigned_taxid"]), "human_host"
                    elif candidate and original_status in {"unclassified", "higher_rank"}:
                        if original_status == "higher_rank" and original_tid:
                            consistent = original_tid in set(taxonomy.lineage(candidate))
                        if consistent:
                            final_tid, source = candidate, "alignment_rescue"
                            accepted += 1
                        else:
                            rejected_lineage += 1
                    info = taxonomy.describe(final_tid)
                    status = ("host_rescued" if source == "human_host" else
                              "alignment_rescued_species" if source == "alignment_rescue" else info["status"])
                    key = (final_tid, status, source)
                    counts[key] += 1
                    metadata[key] = info
                    writer.writerow({
                        "ordinal": n, "r1_read_id": row["r1_read_id"], "r2_read_id": row["r2_read_id"],
                        "assigned_taxid": final_tid, "assigned_rank": info["rank"],
                        "species_taxid": info["species_taxid"], "name": info["name"],
                        "group": info["group"], "status": status, "source": source,
                        "original_taxid": original_tid, "original_status": original_status,
                        "alignment_candidate_taxid": candidate, "lineage_consistent": int(consistent),
                    })
                if n != expected:
                    raise RuntimeError(f"Assignment rows differ from expected: {n} != {expected}")
        finally:
            calls.close()
    return counts, metadata, {"alignment_positive_candidate_fragments": positive,
                              "rescued_fragments": accepted,
                              "rejected_lineage_inconsistent": rejected_lineage}


def materialize_refined(r1, r2, refined_assignments, counts, metadata, out, *, expected, min_fragments,
                        gzip_level=3, max_open=32):
    """Third pass: materialize species meeting the storage threshold as separate FASTQ bins."""
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    species_counts = Counter()
    for (tid, status, source), n in counts.items():
        if status in {"assigned_species", "alignment_rescued_species", "host_rescued"}:
            species_counts[tid] += n
    materialized = {tid for tid, n in species_counts.items() if n >= min_fragments} | ({9606} if 9606 in species_counts else set())
    writers = BufferedWriters(out, max_open=max_open, compresslevel=gzip_level)
    physical = Counter(); tax_rows = Counter(); tax_meta = {}
    stream = fragments(r1, r2)
    completed = False
    try:
        with gzip.open(refined_assignments, "rt", encoding="utf-8", newline="") as inp:
            reader = csv.DictReader(inp, delimiter="\t")
            n = 0
            for n, pair in enumerate(stream, 1):
                row = next(reader, None)
                if row is None or int(row["ordinal"]) != n:
                    raise RuntimeError(f"Refined assignment/FASTQ mismatch at {n}")
                tid = int(row["assigned_taxid"]); status = row["status"]; group = row["group"]
                if status in {"assigned_species", "alignment_rescued_species", "host_rescued"}:
                    base = (f"{group}/{tid}_{safe_name(row['name'])[:100]}" if tid in materialized
                            else f"{group}/species_other")
                elif status == "higher_rank":
                    base = "taxon_unresolved" if group == "unresolved" else f"{group}/taxon_unresolved"
                else:
                    base = "unclassified"
                a, b = pair
                if (row["r1_read_id"] != read_id(a[0]).decode("utf-8", "backslashreplace") or
                        row["r2_read_id"] != read_id(b[0]).decode("utf-8", "backslashreplace")):
                    raise RuntimeError(f"Refined FASTQ read ID differs at ordinal {n}")
                writers.write(fastq_name(base, 1), a)
                writers.write(fastq_name(base, 2), b)
                physical[base] += 1
                key = (tid, status, row["source"], base)
                tax_rows[key] += 1
                tax_meta[key] = row
                if n % 1_000_000 == 0:
                    print(f"refined split {n:,}/{expected:,}", flush=True)
            if n != expected or next(reader, None) is not None:
                raise RuntimeError("Refined assignment/FASTQ total differs")
        completed = True
    finally:
        writers.close(flush=completed); stream.close()
    rows = []
    for key, n in sorted(tax_rows.items(), key=lambda item: (item[0][3], item[0][0])):
        tid, status, source, base = key; row = tax_meta[key]
        rows.append({"assigned_taxid": tid, "assigned_rank": row["assigned_rank"],
                     "species_taxid": row["species_taxid"], "name": row["name"], "group": row["group"],
                     "status": status, "source": source, "bin": base, "fragments": n, "reads": n * 2,
                     "fraction": n / expected, "species_bin_materialized": int(tid in materialized and status in {"assigned_species", "alignment_rescued_species", "host_rescued"})})
    with (out / "counts.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (out / "bin-counts.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(("bin", "fragments", "reads", "fraction"))
        for base, n in sorted(physical.items()): writer.writerow((base, n, n * 2, n / expected))
    return {"materialized_species": len(materialized), "physical_bins": len(physical),
            "fastq_files": len(writers.paths), "species_threshold": min_fragments}


@managed_outputs
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replace", type=replace_bool, default=True)
    p.add_argument("--config")
    p.add_argument("--root")
    p.add_argument("--r1", required=True); p.add_argument("--r2", required=True)
    p.add_argument("--assignments", required=True, help="Existing read-assignments.tsv.gz (v1/v2/v3)")
    p.add_argument("--counts", required=True)
    p.add_argument("--representatives", required=True, help="species-representatives.tsv from metag.sh audit-reference")
    p.add_argument("--out", required=True)
    p.add_argument("--minimap2", help="Default: tools_dir/minimap2/2.28/minimap2, installed by install.sh --qc-tools")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--candidate-min-fragments", type=int, default=50)
    p.add_argument("--max-species", type=int, default=300)
    p.add_argument("--max-reference-gib", type=float, default=20.0)
    p.add_argument("--max-single-reference-gib", type=float, default=4.0)
    p.add_argument("--extra-taxid", action="append", type=int, default=[])
    p.add_argument("--sylph-candidates", help="Optional sylph-candidates.tsv from metag.sh refine --profile-only")
    p.add_argument("--sylph-min-ani", type=float, default=95.0)
    p.add_argument("--include-human", action="store_true",
                   help="Allow taxid 9606/human in generic species refinement; default uses dedicated host validation instead")
    p.add_argument("--rescue-status", action="append", choices=("unclassified", "higher_rank"),
                   default=None, help="Repeat to select eligible statuses; omitted means both")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-aligned-fraction", type=float, default=0.80)
    p.add_argument("--min-identity", type=float, default=0.95)
    p.add_argument("--minimap-index-batch-gib", type=int, default=32,
                   help="Keep the competitive reference in one minimap2 index batch when possible [32]")
    p.add_argument("--materialize", action="store_true", help="Write refined species FASTQs after validation")
    p.add_argument("--materialize-min-fragments", type=int, default=20)
    p.add_argument("--allow-local", action="store_true")
    args = p.parse_args(argv)
    args.rescue_status = args.rescue_status or ["unclassified", "higher_rank"]
    require_compute(args.allow_local)
    if (args.threads < 1 or args.candidate_min_fragments < 1 or args.max_species < 1 or
            args.materialize_min_fragments < 1 or args.minimap_index_batch_gib < 1):
        p.error("thread/species/fragment/index-batch limits must be positive")
    if not 0 <= args.min_mapq <= 255 or not 0 <= args.min_aligned_fraction <= 1 or not 0 <= args.min_identity <= 1:
        p.error("alignment thresholds are outside valid ranges")
    if any(not math.isfinite(x) or x <= 0 for x in (args.max_reference_gib, args.max_single_reference_gib)):
        p.error("reference budgets must be finite and positive")
    if not math.isfinite(args.sylph_min_ani) or not 0 <= args.sylph_min_ani <= 100:
        p.error("--sylph-min-ani must be finite and in 0..100")
    cfg = load_config(args.config, args.root)
    minimap = _executable(args.minimap2 or str(Path(cfg["tools_dir"]) / "minimap2/2.28/minimap2"))
    out = Path(args.out).absolute()
    if not prepare_result(out, replace=args.replace, cfg=cfg,
                          protected=[args.r1, args.r2, args.assignments, args.counts, args.representatives, args.sylph_candidates]):
        return
    check_space(out, float(cfg.get("min_free_gib", 100)))
    reps = _load_representatives(args.representatives, cfg["reference_raw"])
    candidate_counts = _load_candidate_counts(args.counts)
    sylph_extra = _extra_taxids_from_sylph(args.sylph_candidates, args.sylph_min_ani)
    candidate_counts, extra = filter_human_candidates(candidate_counts,
        sorted(set(args.extra_taxid) | set(sylph_extra)), reps, include_human=args.include_human)
    selected, skipped, total_bases = choose_candidates(candidate_counts, reps,
        minimum=args.candidate_min_fragments, maximum=args.max_species,
        max_reference_gib=args.max_reference_gib, max_single_reference_gib=args.max_single_reference_gib,
        extra=extra)
    ref = out / "candidate-reference.fa"
    reference_rows = build_competitive_reference(selected, ref)
    actual_bases = sum(row["bases_written"] for row in reference_rows)
    # minimap2 -I G units are decimal bases; avoid multi-part indexes, whose
    # mapping qualities cannot be interpreted as one competitive reference.
    if actual_bases >= args.minimap_index_batch_gib * 1_000_000_000:
        raise RuntimeError("Candidate reference exceeds one minimap2 index batch; reduce candidates or increase --minimap-index-batch-gib")
    with (out / "candidate-species.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(reference_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(reference_rows)
    with (out / "skipped-candidates.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("species_taxid", "reason", "bases", "kraken_fragments"),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(skipped)
    # Exact expected rows come from assignments; count without reading FASTQ yet.
    expected = 0
    with open_assignment_text(args.assignments) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for expected, row in enumerate(reader, 1):
            if int(row["ordinal"]) != expected:
                raise RuntimeError(f"Non-contiguous assignment ordinal at {expected}")
    if not expected:
        raise RuntimeError("Empty assignment table")
    alignment = align_for_rescue(str(Path(args.r1).resolve(strict=True)), str(Path(args.r2).resolve(strict=True)),
                                 args.assignments, ref, out, expected=expected, minimap2=minimap,
                                 threads=args.threads, rescue_statuses=set(args.rescue_status),
                                 min_mapq=args.min_mapq, min_fraction=args.min_aligned_fraction,
                                 min_identity=args.min_identity, index_batch_gib=args.minimap_index_batch_gib)
    taxonomy = Taxonomy(Path(cfg["metadata_dir"]) / "taxonomy")
    refined = out / "refined-read-assignments.tsv.gz"
    counts, metadata, rescue = apply_rescue(args.assignments, alignment["calls"], taxonomy, expected, refined)
    with (out / "refined-counts.tsv").open("w", newline="", encoding="utf-8") as handle:
        fields = ("assigned_taxid", "name", "group", "status", "source", "fragments")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key, amount in sorted(counts.items()):
            tid, status, source = key
            writer.writerow(dict(assigned_taxid=tid, name=metadata[key]["name"], group=metadata[key]["group"],
                                 status=status, source=source, fragments=amount))
    materialized = None
    if args.materialize:
        materialized = materialize_refined(str(Path(args.r1).resolve(strict=True)), str(Path(args.r2).resolve(strict=True)),
                                           refined, counts, metadata, out / "split", expected=expected,
                                           min_fragments=args.materialize_min_fragments)
    report = {
        "schema": "metag-species-refinement-v1",
        "fragments": expected,
        "candidate_policy": {"minimum_direct_kraken_fragments": args.candidate_min_fragments,
                             "max_species": args.max_species, "max_reference_gib": args.max_reference_gib,
                             "max_single_reference_gib": args.max_single_reference_gib},
        "alignment_policy": {"engine": "minimap2 -x sr", "min_mapq": args.min_mapq,
                             "min_aligned_fraction": args.min_aligned_fraction, "min_identity": args.min_identity,
                             "rescue_statuses": sorted(set(args.rescue_status)),
                             "higher_rank_lineage_guard": True,
                             "mate_conflict_policy": "different positive species => no rescue"},
        "candidate_species": len(selected), "candidate_reference_bases": total_bases,
        "candidate_reference_bases_written": actual_bases,
        "skipped_candidate_count": len(skipped),
        "skipped_candidate_table": str(out / "skipped-candidates.tsv"),
        "candidate_panel_truncated": bool(skipped),
        "skipped_candidates": skipped[:1000],
        "competitive_reference": {"path": str(ref), "sha256": sha256(ref)},
        "alignment": alignment, "rescue": rescue, "materialized": materialized,
        "note": "Original Kraken species calls are preserved. Rescues depend on the selected candidate panel; high MAPQ cannot exclude omitted close relatives. Increased assigned fraction is not validated accuracy."
    }
    atomic_json(out / "species-refinement.json", report)
    finish_result(out)
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
