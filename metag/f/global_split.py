"""Stream original FASTQs into bins using global-classifier ordinal taxids.

The input classification format is exactly one little-endian uint64 taxid per
fragment (one paired fragment or one single read), with no header.  Reads are
never identified by a dictionary of header strings: duplicate names are legal.
The preflight header digest binds their original order to those ordinals.
"""
from __future__ import annotations

from collections import Counter, OrderedDict
from contextlib import ExitStack
import csv
import gzip
import json
from pathlib import Path
import shutil
import struct

from common import DIR_LIB, atomic_json, paired_id, read_id
from fastq_stream import fastq_records, new_header_digest, update_header_digest, validate_pair_ids
from library import file_state
from sample_state import output_state
from split_metacache import Taxonomy, safe_name


ASSIGNMENT_FIELDS = [
    "fragment_index", "fragment_id", "read_id", "mate", "classification_unit",
    "assigned_rank", "assigned_taxid", "species_taxid", "species_name", "group",
    "species_status", "reference_status", "candidate_assembly_accessions",
    "candidate_reference_sequences",
]
REFERENCE_FIELDS = [
    "fragment_index", "fragment_id", "sequence_accession", "reference_file",
    "assembly_accession", "species_taxid", "feature_hits",
]
COUNT_FIELDS = [
    "destination", "group", "species_taxid", "species_name", "assigned_rank",
    "assigned_taxid", "reason", "fragments", "reads",
]
STAGE_SCHEMA = "metag-global-split-stage-v1"
STAGE_MARKER = ".global-split-owner.json"
DISK_CHECK_FRAGMENTS = 100_000


class BufferedWriters:
    """Bound both pending FASTQ bytes and open gzip files independently of reads.

    Small per-bin buffers avoid opening and closing a gzip member for each read
    when hundreds of species are interleaved.  Least-recently-used buffers spill
    before the total budget is exceeded; file handles have a separate LRU bound.
    Eviction creates valid concatenated gzip members when a bin is reopened.
    """

    def __init__(self, root, *, max_open=32, buffer_bytes=64 << 20,
                 file_buffer_bytes=1 << 20):
        if min(max_open, buffer_bytes, file_buffer_bytes) < 1:
            raise ValueError("Writer limits must be positive")
        self.root = Path(root)
        self.max_open = max_open
        self.buffer_bytes = buffer_bytes
        self.file_buffer_bytes = min(file_buffer_bytes, buffer_bytes)
        self.buffers = OrderedDict()
        self.handles = OrderedDict()
        self.pending_bytes = 0
        self.paths = set()

    def _handle(self, key):
        if key in self.handles:
            handle = self.handles.pop(key)
            self.handles[key] = handle
            return handle[0]
        while len(self.handles) >= self.max_open:
            _, (compressed, raw) = self.handles.popitem(last=False)
            try:
                compressed.close()
            finally:
                raw.close()
        path = self.root / key
        path.parent.mkdir(parents=True, exist_ok=True)
        raw = path.open("ab", buffering=256 << 10)
        try:
            compressed = gzip.GzipFile(filename="", fileobj=raw, mode="ab",
                                       compresslevel=1, mtime=0)
        except BaseException:
            raw.close()
            raise
        self.handles[key] = (compressed, raw)
        return compressed

    def _flush(self, key):
        pending = self.buffers.pop(key)
        self.pending_bytes -= len(pending)
        self._handle(key).write(pending)

    def write(self, key, record):
        relative = Path(key)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Bin path must stay within the output directory")
        self.paths.add(key)
        pending = self.buffers.pop(key, bytearray())
        size = sum(map(len, record))
        for line in record:
            pending.extend(line)
        self.pending_bytes += size
        self.buffers[key] = pending
        if len(pending) >= self.file_buffer_bytes:
            self._flush(key)
        while self.pending_bytes > self.buffer_bytes:
            self._flush(next(iter(self.buffers)))

    def close(self, *, flush=True):
        try:
            if flush:
                while self.buffers:
                    self._flush(next(iter(self.buffers)))
            else:
                self.buffers.clear()
                self.pending_bytes = 0
        finally:
            while self.handles:
                _, (compressed, raw) = self.handles.popitem(last=False)
                try:
                    compressed.close()
                finally:
                    raw.close()


def taxids(path, expected_fragments):
    """Read fixed-size blocks; reject truncation or appended assignments."""
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


def prepare_stage(out, input_identity):
    """Restart only a stage bearing this splitter's matching ownership marker."""
    if out.exists() or out.is_symlink():
        raise RuntimeError(f"Output already exists: {out}")
    stage = out.with_name(out.name + ".global.partial")
    owner = {"schema": STAGE_SCHEMA, "input_identity": input_identity,
             "output": str(out.absolute())}
    if stage.is_symlink():
        raise RuntimeError(f"Refusing symlink at split staging path: {stage}")
    if stage.exists():
        try:
            previous = json.loads((stage / STAGE_MARKER).read_text())
        except (OSError, ValueError) as exc:
            raise RuntimeError(f"Unowned split staging path; move it before retrying: {stage}") from exc
        if previous != owner:
            raise RuntimeError(f"Split staging path belongs to different inputs: {stage}")
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    atomic_json(stage / STAGE_MARKER, owner)
    return stage


def check_output_space(path, min_free_bytes):
    """Check the actual destination filesystem, not the scratch filesystem."""
    if min_free_bytes:
        free = shutil.disk_usage(path).free
        if free < min_free_bytes:
            raise RuntimeError(f"Insufficient free split-output space at {path}: "
                               f"{free / 2**30:.2f} GiB free; "
                               f"minimum reserve {min_free_bytes / 2**30:.2f} GiB")


def split(r1, r2, classifications, out, *, expected_fragments, input_identity,
          taxonomy_folder, classifier="Kraken 2 global minimizer classifier",
          before_publish=None, gather=None, expected_header_digest=None,
          min_free_bytes=0, provenance=None, read_reports=True):
    """Publish one complete species split, retaining original FASTQ bytes.

    ``before_publish`` rechecks the root pipeline's input/library identity.
    A failure leaves a marked staging directory that this same run may restart;
    it never publishes split.done in the final output location.
    """
    r1, classifications, out = Path(r1), Path(classifications), Path(out)
    r2 = Path(r2) if r2 is not None else None
    if not isinstance(expected_fragments, int) or expected_fragments < 0:
        raise ValueError("Expected fragment count must be a nonnegative integer")
    if not input_identity:
        raise ValueError("A nonempty input identity is required for split staging")
    if not isinstance(min_free_bytes, int) or min_free_bytes < 0:
        raise ValueError("Minimum free output bytes must be a nonnegative integer")
    classification_state = file_state(classifications)
    if classification_state[0] != 8 * expected_fragments:
        raise RuntimeError("Compact classification size differs from expected fragment count")
    if expected_header_digest is not None and (
            not isinstance(expected_header_digest, str) or len(expected_header_digest) != 64
            or any(c not in "0123456789abcdef" for c in expected_header_digest)):
        raise ValueError("Expected header digest must be a SHA256 hex string")
    tax = Taxonomy(Path(taxonomy_folder))
    out.parent.mkdir(parents=True, exist_ok=True)
    check_output_space(out.parent, min_free_bytes)
    stage = prepare_stage(out, input_identity)
    counts, classification_cache = Counter(), {}
    writers = BufferedWriters(stage)
    digest = new_header_digest()
    a, b = fastq_records(r1), fastq_records(r2) if r2 else None
    n = 0
    scan_complete = False
    try:
        with ExitStack() as stack:
            if read_reports:
                assign_handle = stack.enter_context(gzip.open(
                    stage / "read-assignments.tsv.gz", "wt", compresslevel=1,
                    encoding="utf-8", newline=""))
                assignments = csv.writer(assign_handle, delimiter="\t", lineterminator="\n")
                assignments.writerow(ASSIGNMENT_FIELDS)
                refs_handle = stack.enter_context(gzip.open(
                    stage / "read-reference-hits.tsv.gz", "wt", compresslevel=1,
                    encoding="utf-8", newline=""))
                csv.writer(refs_handle, delimiter="\t", lineterminator="\n").writerow(REFERENCE_FIELDS)
            for tid in taxids(classifications, expected_fragments):
                rec1 = next(a, None)
                rec2 = next(b, None) if b is not None else None
                if rec1 is None or (b is not None and rec2 is None):
                    raise RuntimeError("More compact assignments than FASTQ fragments")
                if rec2:
                    validate_pair_ids(rec1[0], rec2[0])
                update_header_digest(digest, rec1, rec2)
                info = classification_cache.get(tid)
                if info is None:
                    if tid and tid not in tax.parent:
                        raise RuntimeError(f"Compact assignment has unknown taxid: {tid}")
                    rank = tax.rank.get(tid, "none") if tid else "none"
                    species = tax.species(tid)
                    if tid and species and rank != "species":
                        raise RuntimeError(f"Compact assignment is below species rank: {tid}")
                    group = tax.group(tid)
                    species_name = tax.name.get(species, "") if species else ""
                    if rank == "species" and species:
                        base = f"{group}/{species}_{safe_name(species_name or str(species))}"
                        destination, reason, status = "species", "", "assigned_species"
                    else:
                        reason = "unclassified" if tid == 0 else "no_species_in_taxonomic_lineage"
                        base = f"split.bad/{group}/{reason}"
                        destination, status = "split.bad", "higher_rank" if tid else "unclassified"
                    info = (rank, species, group, species_name, base, destination, reason, status)
                    classification_cache[tid] = info
                rank, species, group, species_name, base, destination, reason, status = info
                n += 1
                if read_reports:
                    fragment_id = paired_id(rec1[0], rec2[0]) if rec2 else read_id(rec1[0])
                    for mate, record in ((1, rec1), (2, rec2)) if rec2 else ((0, rec1),):
                        read_name = record[0][1:].split()[0].decode("utf-8", errors="replace")
                        assignments.writerow((
                            n, fragment_id, read_name, mate,
                            "paired_fragment" if r2 else "single_read", rank, tid,
                            species, species_name, group, status, "not_recorded", "[]", "",
                        ))
                if rec2:
                    writers.write(f"{base}/PE.R1.fastq.gz", rec1)
                    writers.write(f"{base}/PE.R2.fastq.gz", rec2)
                else:
                    writers.write(f"{base}/SE.fastq.gz", rec1)
                counts[(destination, group, species, rank, tid, reason)] += 1
                if n % DISK_CHECK_FRAGMENTS == 0:
                    check_output_space(stage, min_free_bytes)
                if n % 1_000_000 == 0:
                    print(f"[split] {n:,}/{expected_fragments:,} fragments; "
                          f"{len(writers.paths):,} FASTQ outputs", flush=True)
            extra1 = next(a, None)
            extra2 = next(b, None) if b is not None else None
            if extra1 is not None or extra2 is not None:
                raise RuntimeError("Fewer compact assignments than FASTQ fragments")
            check_output_space(stage, min_free_bytes)
            scan_complete = True
    finally:
        try:
            # A failed run will restart the stage; do not spend more output
            # space flushing pending read buffers after an error/low-space stop.
            writers.close(flush=scan_complete)
        finally:
            a.close()
            if b is not None:
                b.close()
    if n != expected_fragments:
        raise RuntimeError(f"Split fragment count differs: {n} != {expected_fragments}")
    header_digest = digest.hexdigest()
    if expected_header_digest is not None and header_digest != expected_header_digest:
        raise RuntimeError("FASTQ header/order digest differs from preflight; split was not published")
    rows = []
    for (destination, group, species, rank, tid, reason), fragments in sorted(counts.items()):
        rows.append(dict(destination=destination, group=group, species_taxid=species,
                         species_name=tax.name.get(species, "") if species else "",
                         assigned_rank=rank, assigned_taxid=tid, reason=reason,
                         fragments=fragments, reads=fragments * (2 if r2 else 1)))
    with (stage / "counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COUNT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    reports = {"counts.tsv"}
    if read_reports:
        reports.update({"read-assignments.tsv.gz", "read-reference-hits.tsv.gz"})
    if gather is not None:
        from sample_report import snapshot_references
        snapshot_references(gather, stage / "reference-matches.tsv", DIR_LIB)
        reports.add("reference-matches.tsv")
    if provenance is not None:
        atomic_json(stage / "run-provenance.json", provenance)
        reports.add("run-provenance.json")
    if before_publish is not None:
        before_publish()
    if file_state(classifications) != classification_state:
        raise RuntimeError("Compact assignments changed while splitting; output was not published")
    check_output_space(stage, min_free_bytes)
    if out.exists() or out.is_symlink():
        raise RuntimeError(f"Output appeared while splitting: {out}")
    atomic_json(stage / "split.done", {
        "schema": 2, "input_identity": input_identity,
        "outputs": {name: output_state(stage / name) for name in sorted(writers.paths | reports)},
        "fragments": n, "reads": n * (2 if r2 else 1), "paired": bool(r2), "rows": rows,
        "header_digest": header_digest,
        "reference_status_fragments": {"not_recorded": n},
        "reference_evidence": "not_recorded",
        "sample_reference_evidence": "sourmash_gather" if gather is not None else "not_recorded",
        "reference_evidence_note": "Species-only mode does not collect per-read reference hits; "
                                   "the header-only reference table does not mean that no references matched.",
        "classifier": classifier,
        "read_reports": read_reports,
    })
    # Keep the ownership marker through rename so an interrupted publish remains
    # restartable.  It is a small operational receipt, not a FASTQ/report output.
    stage.rename(out)
    print(f"PASS split: {n:,} fragments -> {out}", flush=True)
    return out
