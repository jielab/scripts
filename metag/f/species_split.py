"""Split original FASTQ records into taxonomic bins with optional human-host rescue.

Key design points
-----------------
1. Kraken taxonomic classification and human-host detection are distinct signals.
2. ``taxid == 0`` (no taxonomic call) is not conflated with a valid higher-rank call.
3. Optional BWA-MEM2 screening against an independent human reference can rescue
   otherwise unresolved vertebrate fragments while preserving the original
   Kraken call in ``read-assignments.tsv.gz``.
4. Both mates always stay together in the final FASTQ bin.

The query context records :func:`split_context` for interpretation. Explicit
--replace FALSE reuses complete outputs; changing policy requires replacement.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack
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

from fastq_stream import fragments, new_digest, read_id, update_digest
from common import atomic_json, check_space, file_state, sha256
from taxonomy import Taxonomy, safe_name
from compressed_output import BufferedWriters, gzip_text
from bin_layout import fastq_name


HOST_INDEX_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
HOST_RESCUE_POLICIES = ("unclassified", "unresolved_vertebrate", "all")
HOST_PAIR_POLICIES = ("any_mate", "both_mates")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")



def taxids(path, expected_fragments):
    """Read compact uint64 assignments and reject truncation/appended records."""
    with Path(path).open("rb", buffering=8 << 20) as handle:
        count = 0
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            if len(block) % 8:
                raise RuntimeError("Truncated compact classification record")
            for (taxid,) in struct.iter_unpack("<Q", block):
                count += 1
                if count > expected_fragments:
                    raise RuntimeError("More compact assignments than expected fragments")
                yield taxid
        if count != expected_fragments:
            raise RuntimeError(f"Compact assignment count differs: {count} != {expected_fragments}")


def _numeric(value, label, *, low=None, high=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label} must be finite")
    if low is not None and value < low:
        raise ValueError(f"{label} must be >= {low}")
    if high is not None and value > high:
        raise ValueError(f"{label} must be <= {high}")
    return value


def _host_reference(cfg):
    value = cfg.get("human_host_reference")
    if not value:
        raise ValueError("human_host_enabled=true requires human_host_reference (indexed FASTA)")
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = Path(cfg.get("root", ".")) / path
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"human_host_reference is not a regular file: {path}")
    return path


def _host_executable(cfg):
    value = str(cfg.get("human_host_executable", "bwa-mem2"))
    candidate = Path(value).expanduser() if os.sep in value else None
    resolved = str(candidate.resolve(strict=True)) if candidate is not None else shutil.which(value)
    if not resolved or not Path(resolved).is_file() or not os.access(resolved, os.X_OK):
        raise ValueError(f"Cannot find executable human aligner: {value!r}")
    return Path(resolved)


def _host_state(path):
    stat = Path(path).stat()
    value = dict(size=stat.st_size,mtime_ns=stat.st_mtime_ns,ctime_ns=stat.st_ctime_ns,inode=stat.st_ino)
    if stat.st_size <= 8 << 20: value["sha256"] = sha256(path)
    return value


def host_context(cfg):
    """Return a small immutable identity record for the optional host screen."""
    enabled = bool(cfg.get("human_host_enabled", False))
    if not enabled:
        return {"enabled": False}
    reference = _host_reference(cfg)
    executable = _host_executable(cfg)
    sidecars = []
    missing = []
    for suffix in HOST_INDEX_SUFFIXES:
        path = Path(str(reference) + suffix)
        if not path.is_file() or not path.stat().st_size:
            missing.append(str(path))
        else:
            stat = path.stat()
            sidecars.append({"path": str(path.resolve()), **_host_state(path)})
    if missing:
        raise ValueError("BWA-MEM2 index is incomplete; run `bwa-mem2 index REF.fa`. Missing: " + ", ".join(missing))
    pair_policy = str(cfg.get("human_host_pair_policy", "any_mate"))
    rescue_policy = str(cfg.get("human_host_rescue_policy", "unresolved_vertebrate"))
    if pair_policy not in HOST_PAIR_POLICIES:
        raise ValueError(f"human_host_pair_policy must be one of {HOST_PAIR_POLICIES}")
    if rescue_policy not in HOST_RESCUE_POLICIES:
        raise ValueError(f"human_host_rescue_policy must be one of {HOST_RESCUE_POLICIES}")
    threads = int(cfg.get("human_host_threads", min(24, int(cfg.get("classify_cores", 24)))))
    if threads < 1 or threads > int(cfg.get("classify_cores", threads)):
        raise ValueError("human_host_threads must be in 1..classify_cores")
    min_mapq = int(cfg.get("human_host_min_mapq", 20))
    if not 0 <= min_mapq <= 255:
        raise ValueError("human_host_min_mapq must be in 0..255")
    min_fraction = float(_numeric(float(cfg.get("human_host_min_aligned_fraction", 0.80)),
                                  "human_host_min_aligned_fraction", low=0, high=1))
    min_identity = float(_numeric(float(cfg.get("human_host_min_identity", 0.95)),
                                  "human_host_min_identity", low=0, high=1))
    regex = cfg.get("human_host_contig_regex")
    if regex is not None:
        re.compile(str(regex))
    excluded = cfg.get("human_host_exclude_contigs", ["chrEBV", "NC_007605.1"])
    if not isinstance(excluded, list) or any(not isinstance(x, str) or not x for x in excluded):
        raise ValueError("human_host_exclude_contigs must be a list of nonempty contig names")
    ref_stat = reference.stat()
    return {
        "enabled": True,
        "engine": "bwa-mem2",
        "executable": str(executable),
        "executable_sha256": sha256(executable),
        "reference": str(reference),
        "reference_state": _host_state(reference),
        "index_files": sidecars,
        "threads": threads,
        "min_mapq": min_mapq,
        "min_aligned_fraction": min_fraction,
        "min_identity": min_identity,
        "pair_policy": pair_policy,
        "rescue_policy": rescue_policy,
        "contig_regex": None if regex is None else str(regex),
        "exclude_contigs": sorted(set(excluded)),
        "meaning": "human-compatible alignment; named viral decoy contigs excluded" if regex is None else
                   "alignment to selected human contigs; species specificity depends on the competitive reference"
    }


def verify_host_context(context):
    if not context.get("enabled"):
        return
    if sha256(context["executable"]) != context["executable_sha256"]:
        raise RuntimeError("Host executable changed during this sample")
    resources = [dict(path=context["reference"], **context["reference_state"]), *context["index_files"]]
    for expected in resources:
        stat = Path(expected["path"]).stat()
        actual = _host_state(expected["path"])
        if any(actual[key] != value for key, value in expected.items() if key != "path"):
            raise RuntimeError(f"Host reference/index changed during this sample: {expected['path']}")


def split_context(cfg):
    """Configuration affecting split output; bind this into query identity.

    ``materialize_species_min_fragments`` controls filesystem materialization,
    not whether a taxon is reported.  Every taxid remains in counts/assignment
    tables; species below the materialization threshold merely share group-level FASTQ files.
    """
    level = int(cfg.get("split_gzip_level", 3))
    if not 1 <= level <= 9:
        raise ValueError("split_gzip_level must be in 1..9")
    minimum = int(cfg.get("materialize_species_min_fragments", 20))
    if minimum < 1:
        raise ValueError("materialize_species_min_fragments must be >=1")
    always = cfg.get("always_materialize_taxids", [9606])
    if not isinstance(always, list) or any(isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in always):
        raise ValueError("always_materialize_taxids must be a list of positive integer taxids")
    return {
        "split_gzip_level": level,
        "write_read_assignments": bool(cfg.get("write_read_assignments", True)),
        "materialize_species_min_fragments": minimum,
        "always_materialize_taxids": sorted(set(always)),
        "human_host": host_context(cfg),
    }


def _stop_process(process):
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


def _cigar_metrics(cigar):
    if cigar == "*":
        return 0, 0
    pos = 0
    aligned_query = query_consumed = 0
    for match in _CIGAR.finditer(cigar):
        if match.start() != pos:
            raise RuntimeError(f"Malformed SAM CIGAR: {cigar!r}")
        length, op = int(match.group(1)), match.group(2)
        if length <= 0:
            raise RuntimeError(f"Malformed SAM CIGAR: {cigar!r}")
        if op in "MIS=X":
            query_consumed += length
        if op in "MI=X":
            aligned_query += length
        pos = match.end()
    if pos != len(cigar):
        raise RuntimeError(f"Malformed SAM CIGAR: {cigar!r}")
    return aligned_query, query_consumed


def _sam_positive(fields, context):
    flag = int(fields[1])
    if flag & 0x4:
        return False, {"mapped": False, "mapq": int(fields[4]), "aligned_fraction": 0.0, "identity": 0.0}
    rname = fields[2]
    # GRCh38 no-alt includes the viral chrEBV decoy. Retain its competition,
    # but do not count alignments to this decoy as human evidence.
    if rname in context.get("exclude_contigs", ["chrEBV", "NC_007605.1"]):
        return False, {"mapped": True, "mapq": int(fields[4]), "aligned_fraction": 0.0, "identity": 0.0}
    if context["contig_regex"] is not None and re.search(context["contig_regex"], rname) is None:
        return False, {"mapped": True, "mapq": int(fields[4]), "aligned_fraction": 0.0, "identity": 0.0}
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
    positive = (mapq >= context["min_mapq"] and
                aligned_fraction >= context["min_aligned_fraction"] and
                identity >= context["min_identity"])
    return positive, {"mapped": True, "mapq": mapq,
                      "aligned_fraction": aligned_fraction, "identity": identity}


def _pair_positive(bits, policy):
    if policy == "any_mate":
        return bool(bits & 0x3)
    if policy == "both_mates":
        return (bits & 0x3) == 0x3
    raise ValueError(f"Unknown pair policy: {policy}")


def screen_human_fastqs(cfg, sample, out, expected_fragments, *, expected_digest=None, limit=0, frozen_context=None):
    """Align a single-end or paired FASTQ stream to the configured host reference.

    Returns a one-byte-per-fragment call file. Bit 0 means R1 passed thresholds,
    bit 1 means R2 passed thresholds. The SAM itself is not retained.
    """
    context = host_context(cfg) if frozen_context is None else frozen_context
    if context != host_context(cfg):
        raise RuntimeError("Host context changed before screening")
    verify_host_context(context)
    if not context["enabled"]:
        return None
    if isinstance(expected_fragments, bool) or not isinstance(expected_fragments, int) or expected_fragments < 1:
        raise ValueError("expected_fragments must be a positive integer")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    calls_path = out / "human-host-calls.bin"
    seen_path = out / ".human-host-seen.bin"
    for path in (calls_path, seen_path):
        if path.exists() or path.is_symlink():
            raise RuntimeError(f"Refusing to overwrite host-screen state: {path}")
    with calls_path.open("w+b") as calls_file, seen_path.open("w+b") as seen_file:
        calls_file.truncate(expected_fragments)
        seen_file.truncate(expected_fragments)
        calls_map = mmap.mmap(calls_file.fileno(), expected_fragments)
        seen_map = mmap.mmap(seen_file.fileno(), expected_fragments)
        result = {}
        errors = []
        stopped = threading.Event()
        command = [context["executable"], "mem"]
        if sample.get("r2") is not None:
            command.append("-p")
        command += ["-t", str(context["threads"]), context["reference"], "-"]
        log_path = out / "human-host.log"
        started = time.monotonic()
        process = producer = None

        def produce():
            digest = new_digest()
            count = 0
            try:
                for a, b in fragments(sample["r1"], sample.get("r2"), limit):
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
                    if sample.get("r2") is not None:
                        if flag & 0x40:
                            bit = 0x1
                        elif flag & 0x80:
                            bit = 0x2
                        else:
                            raise RuntimeError("Paired BWA primary record lacks read1/read2 flag")
                    else:
                        bit = 0x1
                    if seen_map[index] & bit:
                        raise RuntimeError(f"Duplicate primary SAM record for ordinal {ordinal}, bit={bit}")
                    seen_map[index] |= bit
                    primary += 1
                    is_positive, metrics = _sam_positive(fields, context)
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
                expected_seen = 0x3 if sample.get("r2") is not None else 0x1
                missing = sum(1 for i in range(expected_fragments) if seen_map[i] != expected_seen)
                if missing:
                    raise RuntimeError(f"BWA-MEM2 primary SAM coverage incomplete for {missing:,} fragments")
                calls_map.flush()
                calls_file.flush(); os.fsync(calls_file.fileno())
                host_fragments = any_mate = both_mates = r1_pos = r2_pos = 0
                for i in range(expected_fragments):
                    bits = calls_map[i]
                    r1_pos += int(bool(bits & 0x1))
                    r2_pos += int(bool(bits & 0x2))
                    any_mate += int(bool(bits & 0x3))
                    both_mates += int((bits & 0x3) == 0x3)
                    host_fragments += int(_pair_positive(bits, context["pair_policy"]))
                summary = {
                    "schema": "metag-human-host-v1", "status": "complete",
                    "paired": sample.get("r2") is not None,
                    "fragments": expected_fragments, "primary_alignments": primary,
                    "mapped_primary_alignments": mapped,
                    "positive_alignments": positive_alignments,
                    "r1_positive": r1_pos, "r2_positive": r2_pos,
                    "any_mate_positive": any_mate, "both_mates_positive": both_mates,
                    "host_positive_fragments": host_fragments,
                    "host_positive_fraction": host_fragments / expected_fragments,
                    "seconds": time.monotonic() - started,
                    "calls": str(calls_path), "calls_state": file_state(calls_path),
                    "calls_sha256": sha256(calls_path), "context": context,
                    "input_digest": result.get("input_digest")
                }
                verify_host_context(context)
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


def _load_host_calls(path, expected):
    handle = Path(path).open("rb", buffering=8 << 20)
    count = 0
    try:
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            for value in block:
                count += 1
                if count > expected:
                    raise RuntimeError("Host call file has more fragments than expected")
                yield value
        if count != expected:
            raise RuntimeError(f"Host call count differs: {count} != {expected}")
    finally:
        handle.close()


def _rescue_human(tid, info, host_positive, policy):
    if not host_positive or (info["group"] == "human" and info["status"] == "assigned_species"):
        return False
    if policy == "all":
        return True
    if policy == "unclassified":
        return tid == 0
    if policy == "unresolved_vertebrate":
        return tid == 0 or (info["status"] == "higher_rank" and
                            (info.get("human_ancestor", False) or info.get("vertebrate_lineage", False)))
    raise ValueError(f"Unknown human rescue policy: {policy}")


def _base_for(info, tid, policy, materialized_species):
    """Return a disjoint file stem; counts and assignments retain every taxid."""
    group, name = info["group"], info.get("safe_name") or safe_name(info["name"])[:100]
    if info["status"] == "assigned_species":
        if tid in materialized_species or tid in policy["always_materialize_taxids"]:
            return f"{group}/{tid}_{name}"
        return f"{group}/species_other"
    if not tid:
        return "unclassified"
    return "taxon_unresolved" if group == "unresolved" else f"{group}/taxon_unresolved"


def split_originals(cfg, sample, query, snapshot, out, *, host_work=None, frozen_host_context=None, taxonomy=None):
    """Write exact original FASTQ records into evidence-aware physical bins.

    The taxonomy table remains taxid-level even when many low-count taxa share one
    physical FASTQ bin.  This avoids creating thousands of one- or two-fragment gzip files while preserving every assignment.
    """
    split_started = time.monotonic()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    metadata = out / "qc"
    metadata.mkdir()
    reserve = float(cfg.get("min_free_gib", 100))
    check_space(out, reserve)
    policy = split_context(cfg)
    if frozen_host_context is not None and policy["human_host"] != frozen_host_context:
        raise RuntimeError("Host context changed between query and splitting")
    expected = int(query["fragments"])
    classifications = Path(query["taxids"])
    original_state = file_state(classifications)
    if original_state != query["taxids_state"] or original_state["size"] != expected * 8:
        raise RuntimeError("Ordinal classification file changed or is incomplete")
    taxonomy = taxonomy if taxonomy is not None else Taxonomy(snapshot["taxonomy_dir"])

    # A cheap pass over the compact 8-byte/taxon file determines which species
    # meet the file-materialization count threshold.
    # This does not reread the original FASTQ and is tiny relative to index I/O.
    kraken_taxid_counts = Counter(taxids(classifications, expected))
    if sum(kraken_taxid_counts.values()) != expected:
        raise RuntimeError("Pre-counted Kraken assignments do not conserve fragments")
    materialized_species = set()
    species_taxa_total = 0
    species_fragments_total = 0
    for tid, count in kraken_taxid_counts.items():
        info = taxonomy.describe(tid)
        if info["status"] == "assigned_species":
            species_taxa_total += 1
            species_fragments_total += count
            if count >= policy["materialize_species_min_fragments"] or tid in policy["always_materialize_taxids"]:
                materialized_species.add(tid)
    species_fragments_materialized = sum(kraken_taxid_counts[tid] for tid in materialized_species)

    host = None
    if policy["human_host"]["enabled"]:
        if host_work is None:
            import tempfile
            # Direct callers receive an exclusively created scratch directory.
            host_work = Path(tempfile.mkdtemp(prefix=".host-", dir=Path(out).parent))
        host = screen_human_fastqs(cfg, sample, host_work, expected,
                                   expected_digest=query["input_digest"],
                                   limit=int(query["requested_limit"]), frozen_context=policy["human_host"])
        atomic_json(metadata / "human-host-summary.json", host)
        host_iter = _load_host_calls(host["calls"], expected)
    else:
        host_iter = iter(())

    bin_counts = Counter()          # physical FASTQ bins
    taxon_counts = Counter()        # exact final taxid/status/bin rows
    taxon_meta = {}
    kraken_status_counts, kraken_group_counts = Counter(), Counter()
    kraken_group_status = defaultdict(Counter)
    physical_estimates = Counter()
    routes = {}
    for tid, count in kraken_taxid_counts.items():
        info = taxonomy.describe(tid)
        base = _base_for(info, tid, policy, materialized_species)
        separate = info["status"] == "assigned_species" and (tid in materialized_species or tid in policy["always_materialize_taxids"])
        routes[tid] = (info, base, separate)
        physical_estimates[base] += count
    cores = max(1, int(cfg.get("classify_cores", 8)))
    assignment_threads = 2 if cores >= 12 else 1 if cores >= 4 else 0
    hot = [base for base, count in physical_estimates.most_common(2) if count >= 100_000] if cores >= 8 else []
    hot_keys = [fastq_name(base, mate) for base in hot
                for mate in ((1, 2) if sample.get("r2") else (None,))]
    compression_budget = max(1, cores - assignment_threads - 1)
    pigz_threads = min(2, max(1, compression_budget // (len(hot_keys) + 1)))
    compression_threads = min(8, max(1, compression_budget - len(hot_keys) * pigz_threads))
    writers = BufferedWriters(out, max_open=int(cfg.get("max_open_bins", 32)),
                              compresslevel=policy["split_gzip_level"], threads=compression_threads,
                              hot_keys=hot_keys, pigz_threads=pigz_threads)
    digest = new_digest()
    started = time.monotonic()
    n = bases = raw_bytes = 0
    completed = False
    iterator = fragments(sample["r1"], sample.get("r2"), int(query["requested_limit"]))
    host_positive_count = host_rescued_count = host_nonhuman_species_conflicts = 0
    kraken_human_count = 0
    assignment_path = metadata / "read-assignments.tsv.gz"
    assignment_stack = ExitStack()
    assignments = None
    try:
        if policy["write_read_assignments"]:
            assignment_handle = assignment_stack.enter_context(gzip_text(
                assignment_path, compresslevel=policy["split_gzip_level"], threads=assignment_threads))
            assignments = csv.writer(assignment_handle, delimiter="\t", lineterminator="\n")
            assignments.writerow(("ordinal", "r1_read_id", "r2_read_id", "assigned_taxid",
                                  "assigned_rank", "species_taxid", "group", "status", "bin",
                                  "species_bin_materialized",
                                  "kraken_taxid", "kraken_rank", "kraken_group", "kraken_status",
                                  "host_r1", "host_r2", "host_positive", "host_rescued"))
        for kraken_tid in taxids(classifications, expected):
            pair = next(iterator, None)
            if pair is None:
                raise RuntimeError("More ordinal assignments than original FASTQ fragments")
            a, b = pair
            update_digest(digest, a, b)
            kraken_info, original_base, original_separate = routes[kraken_tid]
            if kraken_info["taxid"] != kraken_tid:
                raise RuntimeError("Ordinal taxid was not normalized consistently")
            bits = next(host_iter) if host is not None else 0
            host_positive = (_pair_positive(bits, policy["human_host"].get("pair_policy", "any_mate"))
                             if host is not None else False)
            rescue = (_rescue_human(kraken_tid, kraken_info, host_positive,
                                    policy["human_host"].get("rescue_policy", "unresolved_vertebrate"))
                      if host is not None else False)

            kraken_status_counts[kraken_info["status"]] += 1
            kraken_group_counts[kraken_info["group"]] += 1
            kraken_group_status[kraken_info["group"]][kraken_info["status"]] += 1
            kraken_human_count += int(kraken_info["group"] == "human")
            host_positive_count += int(host_positive)
            host_rescued_count += int(rescue)
            if (host_positive and kraken_info["status"] == "assigned_species" and
                    kraken_info["group"] not in ("human", "unresolved")):
                host_nonhuman_species_conflicts += 1

            if rescue:
                info = {"taxid": 9606, "species_taxid": 9606, "rank": "species",
                        "group": "human", "name": "Homo sapiens", "status": "host_rescued"}
                tid = 9606
                base = "human/9606_Homo_sapiens"
                separately_materialized = True
            else:
                info = kraken_info
                tid = kraken_tid
                base = original_base
                separately_materialized = original_separate

            n += 1
            if assignments is not None:
                assignments.writerow((n, read_id(a[0]).decode("utf-8", "backslashreplace"),
                                      read_id(b[0]).decode("utf-8", "backslashreplace") if b else "",
                                      tid, info["rank"], info["species_taxid"], info["group"],
                                      info["status"], base, int(separately_materialized),
                                      kraken_tid, kraken_info["rank"], kraken_info["group"], kraken_info["status"],
                                      int(bool(bits & 0x1)), int(bool(bits & 0x2)), int(host_positive), int(rescue)))
            if b is None:
                writers.write(fastq_name(base), a)
            else:
                writers.write(fastq_name(base, 1), a)
                writers.write(fastq_name(base, 2), b)
            for record in (a, b) if b is not None else (a,):
                bases += len(record[1].rstrip(b"\r\n"))
                raw_bytes += sum(map(len, record))
            bin_counts[base] += 1
            key = (tid, info["status"], base)
            taxon_counts[key] += 1
            if key not in taxon_meta:
                taxon_meta[key] = {"assigned_taxid": tid, "assigned_rank": info["rank"],
                                   "species_taxid": info["species_taxid"], "name": info["name"],
                                   "group": info["group"], "status": info["status"], "bin": base,
                                   "species_bin_materialized": int(separately_materialized)}
            if n % 100_000 == 0:
                check_space(out, reserve)
            if n % 1_000_000 == 0:
                print(f"[{sample['sample_id']}] split {n:,}/{expected:,}", flush=True)
        if next(iterator, None) is not None:
            raise RuntimeError("Fewer ordinal assignments than original FASTQ fragments")
        if host is not None:
            try:
                next(host_iter)
                raise RuntimeError("More host calls than expected fragments")
            except StopIteration:
                pass
        completed = True
    finally:
        try:
            writers.close(flush=completed)
        finally:
            iterator.close()
            assignment_stack.close()

    if (n != expected or digest.hexdigest() != query["input_digest"] or
            bases != query["total_bases"] or raw_bytes != query["fastq_uncompressed_bytes"]):
        raise RuntimeError("Original FASTQ bytes/order/count differ from classification pass")
    if file_state(classifications) != original_state:
        raise RuntimeError("Ordinal classification file changed while splitting")

    paired_multiplier = 2 if sample.get("r2") else 1
    rows = []
    for key in sorted(taxon_counts, key=lambda k: (taxon_meta[k]["bin"], k[0], k[1])):
        row = dict(taxon_meta[key])
        value = taxon_counts[key]
        row.update(fragments=value, reads=value * paired_multiplier,
                   fraction=value / expected)
        rows.append(row)
    with (metadata / "counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

    bin_taxids, bin_species = defaultdict(set), defaultdict(set)
    for row in rows:
        bin_taxids[row["bin"]].add(row["assigned_taxid"])
        if row["species_taxid"]:
            bin_species[row["bin"]].add(row["species_taxid"])
    bin_rows = [{"bin": base, "fragments": count, "reads": count * paired_multiplier,
                 "fraction": count / expected,
                 "assigned_taxids": ",".join(str(tid) for tid in sorted(bin_taxids[base])),
                 "species_taxids": ",".join(str(tid) for tid in sorted(bin_species[base]))}
                for base, count in sorted(bin_counts.items())]
    with (metadata / "bin-counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("bin", "fragments", "reads", "fraction", "assigned_taxids", "species_taxids"),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(bin_rows)

    final_status_counts, final_group_counts = Counter(), Counter()
    for row in rows:
        final_status_counts[row["status"]] += row["fragments"]
        final_group_counts[row["group"]] += row["fragments"]
    if sum(bin_counts.values()) != expected or sum(row["reads"] for row in rows) != query["reads"]:
        raise RuntimeError("Output count conservation check failed")

    species_taxa_materialized = len(materialized_species)
    species_other_fragments = species_fragments_total - species_fragments_materialized
    result = {
        "schema": "metag-hpc-split-v4", "sample_id": sample["sample_id"],
        "fragments": n, "reads": query["reads"], "paired": query["paired"],
        "total_bases": bases, "input_digest": digest.hexdigest(),
        "validated_scope": query["validated_scope"], "requested_limit": query["requested_limit"],
        "classification_seconds": query["seconds"], "split_seconds": time.monotonic() - split_started,
        "split_preparation_seconds": started - split_started,
        "groups": dict(final_group_counts), "statuses": dict(final_status_counts),
        "kraken_groups": dict(kraken_group_counts), "kraken_statuses": dict(kraken_status_counts),
        "kraken_group_status": {group: dict(values) for group, values in sorted(kraken_group_status.items())},
        "kraken_unclassified_fragments": kraken_status_counts.get("unclassified", 0),
        "kraken_unclassified_fraction": kraken_status_counts.get("unclassified", 0) / expected,
        "kraken_higher_rank_fragments": kraken_status_counts.get("higher_rank", 0),
        "kraken_higher_rank_fraction": kraken_status_counts.get("higher_rank", 0) / expected,
        "kraken_species_fragments": kraken_status_counts.get("assigned_species", 0),
        "kraken_species_fraction": kraken_status_counts.get("assigned_species", 0) / expected,
        "kraken_classified_fragments": expected - kraken_status_counts.get("unclassified", 0),
        "kraken_classified_fraction": 1 - kraken_status_counts.get("unclassified", 0) / expected,
        "kraken_human_fragments": kraken_human_count,
        "host_positive_fragments": host_positive_count,
        "host_rescued_fragments": host_rescued_count,
        "host_positive_nonhuman_species_conflicts": host_nonhuman_species_conflicts,
        "final_human_fragments": final_group_counts.get("human", 0),
        "species_taxa_total": species_taxa_total,
        "species_taxa_materialized": species_taxa_materialized,
        "species_fragments_materialized": species_fragments_materialized,
        "species_fragments_materialized_fraction_of_species_calls": (
            species_fragments_materialized / species_fragments_total if species_fragments_total else None),
        "species_other_fragments": species_other_fragments,
        "physical_bins": len(bin_counts),
        "fastq_files": len(writers.paths), "count_conservation": True,
        "pair_policy": "both original mates share one final fragment/bin; host screen may evaluate mates independently",
        "taxonomy_policy": "exact taxids remain in tables; higher-rank and true unclassified calls are separated; species below the materialization threshold may share physical FASTQ bins",
        "output_bins": {
            "TAXID_name": "one species with at least the materialization threshold; original mates stay together",
            "species_other": "species below the file-materialization threshold within this group; exact species retained in tables",
            "taxon_unresolved": "valid higher-rank calls without a species; root means no unique taxonomic group",
            "unclassified": "taxid 0, no retained taxonomic assignment; not a confidence grade"
        },
        "compression": {"gzip_threads": compression_threads, "pigz_large_streams": len(writers.hot_keys),
                        "pigz_threads_per_stream": pigz_threads, "chunk_bytes": writers.file_buffer_bytes},
        "batch_samples": query.get("batch_samples", 1),
        "batch_index_load_seconds": query.get("batch_index_load_seconds"),
        "batch_taxonomy_load_seconds": query.get("batch_taxonomy_load_seconds"),
        "project_preparation": query.get("project_preparation", {}),
        "batch_query_seconds_excluding_index_load": query.get("batch_query_seconds_excluding_index_load"),
        "batch_classification_seconds": query.get("batch_classification_seconds", query["seconds"]),
        "classification_timing": query.get("classification_timing", "sample wall time"),
        "split_policy": policy,
        "reference_hits": "Kraken taxonomic calls plus optional independent BWA-MEM2 human-host alignment"
    }
    atomic_json(metadata / "summary.json", result)
    check_space(out, reserve)
    return result

