"""Derive theoretical ALT candidates only when an annotation command needs them.

The input run must already have passed CLI load_run/verify_shard. Original
shards are read-only. Committed eager tables remain usable without rebuilding;
deferred tables are reconstructed with bounded memory from matching site and
allele groups. Raw formatted AC/AN/AF strings are copied, never recomputed.

Derived outputs use the workflow's explicitly labelled size/mtime/edge hash
fingerprints, not a claimed full-file checksum. Each shard is committed
atomically and verified on reuse; interrupted shards can be regenerated while
completed ones are reused. A completed but modified cache is never overwritten.
"""
from __future__ import annotations

from contextlib import ExitStack
import csv
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import tempfile

from core import dump_json, fingerprint, open_text, stable_hash, tsv_writer
from scan import BASES, CANDIDATE_FIELDS, relation

SCHEMA = 1
CANDIDATES = "candidates.tsv.gz"
SITES = "multiallelic_sites.tsv.gz"
ALLELES = "multiallelic_alleles.tsv.gz"
SITE_REQUIRED = {"chrom", "pos", "ref", "context", "cpg_site", "n_alt", "primary_alt", "primary_tie"}
ALLELE_REQUIRED = {"chrom", "pos", "ref", "alt", "status", "is_primary", "ac", "an", "af"}
PRESENT_STATUSES = {"observed_pass", "present_failed_filter", "present_missing_counts", "present_below_threshold"}


def _rows(path, required):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not required.issubset(fields):
            raise ValueError(f"Invalid candidate source schema: {path}; required {sorted(required)}")
        for row in reader:
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f"Malformed candidate source row at {path}:{reader.line_num}")
            yield row


def _coordinate(row):
    try:
        pos = int(row["pos"])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid candidate source position: {row.get('pos')}") from exc
    if pos < 1 or not row["chrom"] or len(row["ref"]) != 1 or row["ref"] not in BASES:
        raise ValueError("Candidate sources require positive positions and A/C/G/T REF")
    return row["chrom"], pos, row["ref"]


def _write_candidates(site_path, allele_path, output_path, expected_rows):
    """Join the two already synchronized shard tables, never a genome-size dict."""
    with ExitStack() as stack:
        sites = _rows(site_path, SITE_REQUIRED)
        alleles = _rows(allele_path, ALLELE_REQUIRED)
        stack.callback(sites.close)
        stack.callback(alleles.close)
        groups = itertools.groupby(alleles, key=_coordinate)
        handle, writer = tsv_writer(output_path, CANDIDATE_FIELDS)
        stack.enter_context(handle)
        rows_written = sites_written = 0
        previous = None
        completed_chroms = set()
        for site in sites:
            key = _coordinate(site)
            chrom, pos, ref = key
            if previous:
                if chrom == previous[0] and pos <= previous[1]:
                    raise ValueError(f"Duplicate or unsorted source site: {chrom}:{pos}")
                if chrom != previous[0]:
                    completed_chroms.add(previous[0])
                    if chrom in completed_chroms:
                        raise ValueError(f"Non-contiguous source chromosome: {chrom}")
            previous = key
            try:
                allele_key, group = next(groups)
            except StopIteration as exc:
                raise ValueError(f"Missing allele group for candidate site {key}") from exc
            if key != allele_key:
                raise ValueError(f"Site/allele candidate sources disagree: {key} versus {allele_key}")
            entries = {}
            for allele in group:
                alt = allele["alt"]
                if len(alt) != 1 or alt not in BASES or alt == ref or alt in entries:
                    raise ValueError(f"Invalid or duplicate source ALT at {key}: {alt}")
                if allele["status"] not in PRESENT_STATUSES:
                    raise ValueError(f"Unknown source ALT status at {key}: {allele['status']}")
                if allele["is_primary"] != str(int(alt == site["primary_alt"])):
                    raise ValueError(f"Source primary ALT labels disagree at {key}")
                entries[alt] = allele
            observed = {alt for alt, entry in entries.items() if entry["status"] == "observed_pass"}
            if site["n_alt"] not in ("2", "3") or len(observed) != int(site["n_alt"]):
                raise ValueError(f"Source ALT occupancy disagrees at {key}")
            if site["primary_alt"] not in observed or site["primary_tie"] not in ("0", "1"):
                raise ValueError(f"Invalid primary ALT/tie at {key}")
            if len(site["context"]) != 3 or site["context"][1] != ref or site["cpg_site"] not in ("0", "1"):
                raise ValueError(f"Invalid source context/CpG flag at {key}")
            for alt in BASES:
                if alt == ref:
                    continue
                entry = entries.get(alt)
                transition = int(relation(ref, alt) == "transition")
                writer.writerow(dict(chrom=site["chrom"], pos=site["pos"], ref=ref, alt=alt,
                    status=entry["status"] if entry else "not_listed_in_release",
                    is_primary=int(alt == site["primary_alt"]), primary_tie=site["primary_tie"],
                    n_alt=site["n_alt"], ac=entry["ac"] if entry else ".",
                    an=entry["an"] if entry else ".", af=entry["af"] if entry else ".",
                    context=site["context"], primary_alt=site["primary_alt"],
                    transition=transition, cpg_transition=int(site["cpg_site"] == "1" and transition)))
                rows_written += 1
            sites_written += 1
        if next(groups, None) is not None:
            raise ValueError("Unmatched trailing allele groups in candidate sources")
        if rows_written != expected_rows:
            raise ValueError(f"Candidate row count differs from scan: expected {expected_rows}, wrote {rows_written}")
    return dict(rows=rows_written, sites=sites_written)


def _verified_source(shard, marker, name):
    expected = marker.get("outputs", {}).get(name)
    path = shard / name
    if not isinstance(expected, dict) or not path.is_file():
        raise ValueError(f"Missing committed candidate source file: {path}")
    actual = fingerprint(path)
    if any(actual[key] != expected.get(key) for key in ("size", "edge_sha256")):
        raise ValueError(f"Committed candidate source changed: {path}")
    return actual


def _expected_rows(scan):
    expected = scan.get("expected_candidate_rows")
    n_multi = scan.get("qc", {}).get("multi_alt_sites")
    if expected is None and type(n_multi) is int:
        expected = 3*n_multi  # Old eager scans recorded this QC count.
    if expected is not None and (type(expected) is not int or expected < 0 or expected % 3):
        raise ValueError("Invalid expected_candidate_rows in scan metadata")
    if type(n_multi) is int and expected is not None and expected != 3*n_multi:
        raise ValueError("Candidate row expectation disagrees with multi-ALT scan count")
    return expected


def _verify_derived(target, signature, identity, expected_rows):
    marker = target / "COMPLETE.json"
    if not marker.is_file():
        raise ValueError(f"Unrecognized/incomplete derived candidate cache preserved at {target}; use a new cache directory")
    saved = json.loads(marker.read_text())
    if (saved.get("complete") is not True or saved.get("signature") != signature
            or saved.get("identity") != identity or saved.get("rows") != expected_rows):
        raise ValueError(f"Derived candidate cache identity/count changed: {target}")
    path = target / CANDIDATES
    if not path.is_file():
        raise ValueError(f"Derived candidate cache output is missing: {path}")
    actual, expected = fingerprint(path), saved.get("output", {})
    if any(actual[key] != expected.get(key) for key in ("size", "mtime_ns", "edge_sha256")):
        raise ValueError(f"Derived candidate cache output changed: {path}; preserve it and use a new cache directory")
    return path


def candidate_files_for_run(run, out_dir):
    """Return reusable candidate paths for a previously verified local scan run.

    ``out_dir`` is a separate derived-cache directory, e.g. OUT/candidate_cache.
    It is not created when all shards already contain committed eager tables.
    No original shard, completed cache, or unknown pre-existing output is removed.
    """
    if run.get("status") not in ("scanned", "complete") or not run.get("signature"):
        raise ValueError("candidate_files_for_run requires a verified scanned/complete run")
    shards = [Path(p).expanduser().resolve() for p in run["shards"]]
    if len(set(shards)) != len(shards):
        raise ValueError("Duplicate shards in candidate source run")
    cache_root = Path(out_dir).expanduser().resolve()
    if any(cache_root == shard or cache_root.is_relative_to(shard) for shard in shards):
        raise ValueError("Derived candidate cache must be outside original scan shards")
    generator_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(files=[], source_shards=[], expected_candidate_rows=0,
                  reused_eager_shards=0, generated_shards=0, reused_derived_shards=0,
                  cache_dir=str(cache_root), schema=SCHEMA,
                  validation="metadata/edge fingerprints; not full-file cryptographic checksums")
    all_expected_known = True
    for shard in shards:
        marker_path = shard / "COMPLETE.json"
        if not marker_path.is_file():
            raise ValueError(f"Missing original shard completion marker: {marker_path}")
        marker = json.loads(marker_path.read_text())
        if marker.get("signature") != run["signature"]:
            raise ValueError(f"Original shard/run identity mismatch at {shard}")
        scan_fingerprint = _verified_source(shard, marker, "scan.json")
        scan = json.loads((shard / "scan.json").read_text())
        expected_rows = _expected_rows(scan)
        if expected_rows is None:
            all_expected_known = False
        else:
            result["expected_candidate_rows"] += expected_rows
        if CANDIDATES in marker.get("outputs", {}):
            if scan.get("candidate_storage", "eager") != "eager":
                raise ValueError(f"Candidate storage metadata contradicts committed eager table at {shard}")
            _verified_source(shard, marker, CANDIDATES)
            candidate_path = shard / CANDIDATES
            mode = "reused_eager"
            result["reused_eager_shards"] += 1
            cache_manifest = None
        else:
            if scan.get("candidate_storage") != "deferred":
                raise ValueError(f"No committed eager candidate table and no explicit deferred mode at {shard}")
            if expected_rows is None:
                raise ValueError(f"Deferred scan lacks expected candidate count at {shard}")
            sources = {name: _verified_source(shard, marker, name) for name in (SITES, ALLELES)}
            sources["scan.json"] = scan_fingerprint
            identity = dict(schema=SCHEMA, columns=CANDIDATE_FIELDS, generator_sha256=generator_sha256,
                            source_run_signature=run["signature"], source_shard=str(shard),
                            sources=sources, expected_rows=expected_rows)
            signature = stable_hash(identity)
            cache_root.mkdir(parents=True, exist_ok=True)
            target = cache_root / signature
            with (cache_root / (signature + ".lock")).open("a") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                if target.exists():
                    candidate_path = _verify_derived(target, signature, identity, expected_rows)
                    mode = "reused_derived"
                    result["reused_derived_shards"] += 1
                else:
                    # A new uniquely owned staging directory means hard-kill
                    # leftovers are preserved, not adopted or silently deleted.
                    with tempfile.TemporaryDirectory(prefix=signature + ".partial-", dir=cache_root) as temporary:
                        stage = Path(temporary)
                        stats = _write_candidates(shard / SITES, shard / ALLELES, stage / CANDIDATES, expected_rows)
                        if any(fingerprint(shard / name) != before for name, before in sources.items()):
                            raise ValueError(f"Original candidate sources changed during reconstruction: {shard}")
                        output = fingerprint(stage / CANDIDATES)
                        output["path"] = str(target / CANDIDATES)
                        dump_json(stage / "COMPLETE.json", dict(signature=signature, identity=identity,
                                  output=output, rows=stats["rows"], sites=stats["sites"], complete=True))
                        os.replace(stage, target)
                    candidate_path = target / CANDIDATES
                    mode = "generated"
                    result["generated_shards"] += 1
            cache_manifest = str(target / "COMPLETE.json")
        result["files"].append(str(candidate_path))
        result["source_shards"].append(dict(shard=str(shard), mode=mode,
                    candidate_file=str(candidate_path), expected_rows=expected_rows, cache_manifest=cache_manifest))
    if not all_expected_known:
        result["expected_candidate_rows"] = None
    return result
