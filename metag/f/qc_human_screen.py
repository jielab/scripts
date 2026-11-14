"""Independent human QC alignment, separate from production split semantics.

This module preserves the enhanced candidate-bin streaming and alignment
inspection API while species_split.py remains the production implementation.
The stable host-reference, pair-policy and FASTQ helpers are reused unchanged.
"""
from __future__ import annotations

import json
import mmap
import os
from pathlib import Path
import re
import subprocess
import threading
import time

from common import atomic_json, file_state, sha256
from fastq_stream import fragments, new_digest, update_digest
from species_split import _cigar_metrics, _pair_positive, _stop_process, host_context

__all__ = ["screen_human_fastqs", "host_context", "_pair_positive", "_sam_positive"]


def _sam_positive(fields, context):
    flag = int(fields[1])
    if flag & 0x4:
        return False, {"mapped": False, "mapq": int(fields[4]), "aligned_fraction": 0.0, "identity": 0.0}
    rname = fields[2]
    # GRCh38 no-alt includes the viral chrEBV decoy. Retain its competition,
    # but do not count alignments to this decoy as human evidence.
    excluded = (rname in context.get("exclude_contigs", ["chrEBV", "NC_007605.1"]) or
                (context["contig_regex"] is not None and re.search(context["contig_regex"], rname) is None))
    mapq = int(fields[4])
    aligned, consumed = _cigar_metrics(fields[5])
    sequence = fields[9]
    read_length = len(sequence) if sequence != "*" else consumed
    if read_length <= 0:
        raise RuntimeError("Primary SAM alignment has no query length")
    aligned_fraction = aligned / read_length
    nm = None
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            nm = int(field[5:])
            break
    if nm is None:
        raise RuntimeError("Primary BWA-MEM2 SAM record lacks NM:i tag")
    identity = max(0.0, 1.0 - nm / max(1, aligned))
    positive = (not excluded and mapq >= context["min_mapq"] and
                aligned_fraction >= context["min_aligned_fraction"] and
                identity >= context["min_identity"])
    return positive, {"mapped": True, "mapq": mapq,
                      "aligned_fraction": aligned_fraction, "identity": identity}


def screen_human_fastqs(cfg, sample, out, expected_fragments, *, expected_digest=None, limit=0,
                       input_sources=None, evidence_database=None):
    """Align a single-end or paired FASTQ stream to the configured host reference.

    Returns a one-byte-per-fragment call file. Bit 0 means R1 passed thresholds,
    bit 1 means R2 passed thresholds. The SAM itself is not retained. Optional
    input_sources streams multiple candidate bins through ONE aligner process;
    all sources must use the same single/paired mode.
    """
    context = host_context(cfg)
    if not context["enabled"]:
        return None
    if isinstance(expected_fragments, bool) or not isinstance(expected_fragments, int) or expected_fragments < 1:
        raise ValueError("expected_fragments must be a positive integer")
    sources = list(input_sources) if input_sources is not None else [sample]
    if not sources or any(not source.get("r1") for source in sources):
        raise ValueError("Host-screen input sources must contain R1 paths")
    paired = sample.get("r2") is not None
    if any(bool(source.get("r2")) != paired for source in sources):
        raise ValueError("Host-screen sources cannot mix single-end and paired-end inputs")
    if limit and len(sources) != 1:
        raise ValueError("A fragment limit is only supported for a single host-screen input")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    calls_path = out / "human-host-calls.bin"
    seen_path = out / ".human-host-seen.bin"
    for path in (calls_path, seen_path):
        if path.exists() or path.is_symlink():
            raise RuntimeError(f"Refusing to overwrite host-screen state: {path}")
    evidence = None
    evidence_rows = []
    if evidence_database is not None:
        import sqlite3
        evidence_path = Path(evidence_database)
        if evidence_path.exists():
            raise RuntimeError(f"Refusing to overwrite alignment evidence: {evidence_path}")
        evidence = sqlite3.connect(str(evidence_path))
        evidence.execute("PRAGMA journal_mode=OFF")
        evidence.execute("PRAGMA synchronous=OFF")
        evidence.execute("CREATE TABLE alignments (ordinal INTEGER, mate INTEGER, flag INTEGER, contig TEXT, "
                         "position INTEGER, cigar TEXT, mapq INTEGER, aligned_fraction REAL, identity REAL, "
                         "passed INTEGER, PRIMARY KEY (ordinal, mate))")
    with calls_path.open("w+b") as calls_file, seen_path.open("w+b") as seen_file:
        calls_file.truncate(expected_fragments)
        seen_file.truncate(expected_fragments)
        calls_map = mmap.mmap(calls_file.fileno(), expected_fragments)
        seen_map = mmap.mmap(seen_file.fileno(), expected_fragments)
        result = {}
        errors = []
        stopped = threading.Event()
        command = [context["executable"], "mem"]
        if paired:
            command.append("-p")
        command += ["-t", str(context["threads"]), context["reference"], "-"]
        log_path = out / "human-host.log"
        started = time.monotonic()
        process = producer = None

        def input_fragments():
            for source in sources:
                yield from fragments(source["r1"], source.get("r2"), limit)

        def produce():
            digest = new_digest()
            count = 0
            try:
                for a, b in input_fragments():
                    if stopped.is_set():
                        raise RuntimeError("Human host FASTQ stream interrupted")
                    count += 1
                    if count > expected_fragments:
                        raise RuntimeError("Host-screen input has more fragments than Kraken query")
                    update_digest(digest, a, b)
                    records = (a, b) if b is not None else (a,)
                    qname = str(count).encode("ascii")
                    for record in records:
                        process.stdin.write(b"@" + qname + b"\n" + record[1].rstrip(b"\r\n") +
                                            b"\n+\n" + record[3].rstrip(b"\r\n") + b"\n")
                    if count % 1_000_000 == 0:
                        print(f"[{sample.get('sample_id', sample.get('sample', 'sample'))}] human screen streamed {count:,}/{expected_fragments:,}", flush=True)
                if count != expected_fragments:
                    raise RuntimeError(f"Host-screen FASTQ count differs: {count} != {expected_fragments}")
                if expected_digest is not None and digest.hexdigest() != expected_digest:
                    raise RuntimeError("Host-screen FASTQ digest differs from Kraken query input")
                result.update(fragments=count, input_digest=digest.hexdigest())
                process.stdin.close()
            except BaseException as exc:
                errors.append(exc)
                stopped.set()
                _stop_process(process)
                try:
                    process.stdin.close()
                except (OSError, ValueError):
                    pass

        try:
            with log_path.open("w", encoding="utf-8") as log:
                log.write("COMMAND " + json.dumps(command) + "\n")
                log.write("HOST_CONTEXT " + json.dumps(context, sort_keys=True) + "\n")
                log.flush()
                process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                           stderr=log, start_new_session=True, bufsize=1 << 20)
                producer = threading.Thread(target=produce, name="metag-human-producer", daemon=True)
                producer.start()
                primary = mapped = positive_alignments = 0
                for raw in process.stdout:
                    if raw.startswith(b"@"):
                        continue
                    fields_b = raw.rstrip(b"\r\n").split(b"\t")
                    if len(fields_b) < 11:
                        raise RuntimeError(f"Malformed SAM row: {raw[:160]!r}")
                    fields = [x.decode("utf-8", "strict") for x in fields_b]
                    try:
                        ordinal = int(fields[0])
                        flag = int(fields[1])
                    except ValueError as exc:
                        raise RuntimeError(f"Invalid numeric SAM qname/flag: {fields[:2]}") from exc
                    if not 1 <= ordinal <= expected_fragments:
                        raise RuntimeError(f"SAM ordinal outside expected range: {ordinal}")
                    if flag & (0x100 | 0x800):
                        continue
                    index = ordinal - 1
                    if paired:
                        mate_flags = flag & (0x40 | 0x80)
                        if not flag & 0x1 or mate_flags not in (0x40, 0x80):
                            raise RuntimeError("Paired BWA primary record has invalid paired/read1/read2 flags")
                        if mate_flags == 0x40:
                            bit = 0x1
                        else:
                            bit = 0x2
                    else:
                        bit = 0x1
                    if seen_map[index] & bit:
                        raise RuntimeError(f"Duplicate primary SAM record for ordinal {ordinal}, bit={bit}")
                    seen_map[index] |= bit
                    if paired and flag & 0x2:
                        seen_map[index] |= bit << 2
                    primary += 1
                    is_positive, metrics = _sam_positive(fields, context)
                    if evidence is not None:
                        evidence_rows.append((ordinal, bit, flag, fields[2], int(fields[3]), fields[5],
                                              metrics["mapq"], metrics["aligned_fraction"], metrics["identity"], int(is_positive)))
                        if len(evidence_rows) >= 4096:
                            evidence.executemany("INSERT INTO alignments VALUES (?,?,?,?,?,?,?,?,?,?)", evidence_rows)
                            evidence_rows.clear()
                    mapped += int(metrics["mapped"])
                    if is_positive:
                        calls_map[index] |= bit
                        positive_alignments += 1
                code = process.wait()
                producer.join(timeout=15)
                if producer.is_alive():
                    raise RuntimeError("Human host FASTQ producer did not stop")
                if errors:
                    raise RuntimeError(f"Human host FASTQ producer failed: {errors[0]}") from errors[0]
                if code:
                    raise RuntimeError(f"BWA-MEM2 exited {code}; see {log_path}")
                expected_seen = 0x3 if paired else 0x1
                missing = sum(1 for i in range(expected_fragments) if seen_map[i] & 0x3 != expected_seen)
                if missing:
                    raise RuntimeError(f"BWA-MEM2 primary SAM coverage incomplete for {missing:,} fragments")
                contig_counts = None
                if evidence is not None:
                    if evidence_rows:
                        evidence.executemany("INSERT INTO alignments VALUES (?,?,?,?,?,?,?,?,?,?)", evidence_rows)
                        evidence_rows.clear()
                    evidence.commit()
                    contig_counts = dict(evidence.execute("SELECT contig, COUNT(*) FROM alignments WHERE passed=1 GROUP BY contig ORDER BY contig"))
                calls_map.flush()
                calls_file.flush(); os.fsync(calls_file.fileno())
                host_fragments = any_mate = both_mates = concordant_pairs = r1_pos = r2_pos = 0
                for i in range(expected_fragments):
                    bits = calls_map[i]
                    r1_pos += int(bool(bits & 0x1))
                    r2_pos += int(bool(bits & 0x2))
                    any_mate += int(bool(bits & 0x3))
                    both_mates += int((bits & 0x3) == 0x3)
                    concordant_pairs += int((bits & 0x3) == 0x3 and (seen_map[i] & 0xc) == 0xc)
                    host_fragments += int(_pair_positive(bits, context["pair_policy"]))
                summary = {
                    "schema": "metag-human-host-v1", "status": "complete",
                    "paired": paired,
                    "fragments": expected_fragments, "primary_alignments": primary,
                    "mapped_primary_alignments": mapped,
                    "positive_alignments": positive_alignments,
                    "r1_positive": r1_pos, "r2_positive": r2_pos,
                    "any_mate_positive": any_mate, "both_mates_positive": both_mates,
                    "exactly_one_mate_positive": any_mate - both_mates if paired else None,
                    "neither_mate_positive": expected_fragments - any_mate,
                    "both_mates_concordant_positive": concordant_pairs if paired else None,
                    "host_positive_fragments": host_fragments,
                    "host_positive_fraction": host_fragments / expected_fragments,
                    "seconds": time.monotonic() - started,
                    "calls": str(calls_path), "calls_state": file_state(calls_path),
                    "calls_sha256": sha256(calls_path), "context": context,
                    "input_digest": result.get("input_digest")
                }
                if contig_counts is not None:
                    summary["positive_alignment_contig_counts"] = contig_counts
                    summary["positive_alignment_distinct_contigs"] = len(contig_counts)
                atomic_json(out / "human-host-summary.json", summary)
                seen_path.unlink(missing_ok=True)
                return summary
        finally:
            stopped.set()
            _stop_process(process)
            if producer is not None:
                producer.join(timeout=10)
            if process is not None:
                if process.stdout:
                    process.stdout.close()
                try:
                    process.stdin.close()
                except (OSError, ValueError, AttributeError):
                    pass
            calls_map.close()
            seen_map.close()
            if evidence is not None:
                evidence.close()
