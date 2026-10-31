#!/usr/bin/env python3
from __future__ import annotations
import csv
import hashlib
import io
import json
import re
import sys
from pathlib import Path
from common import DIR_REF, DIR_META, DIR_LIB, REF_GROUPS, reference_files, atomic_json, assert_reference_policy
from local_metadata import TAXONOMY_FILES, prepare_taxonomy, summary_candidates, summary_path, taxonomy_path

SUMMARY_URLS = {
    "viral": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt",
    "archaea": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt",
    "bacteria": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
    "fungi": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt",
    "protozoa": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/assembly_summary.txt",
    "plant": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt",
    "invertebrate": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/assembly_summary.txt",
    "vertebrate_mammalian": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/assembly_summary.txt",
    "vertebrate_other": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/assembly_summary.txt",
}

ACC_RE = re.compile(r"^(GC[AF]_\d+\.\d+)_")

def accession_from_filename(path: Path):
    m = ACC_RE.match(path.name)
    return m.group(1) if m else None


def read_metadata_source(source):
    marker = source.with_name(source.name + ".aria2")
    if marker.exists():
        raise RuntimeError(f"Metadata download is incomplete: {source}; finish aria2c before building")
    raw = source.read_bytes()
    if marker.exists():
        raise RuntimeError(f"Metadata changed during download: {source}; finish aria2c before building")
    return raw

def check_metadata(include_taxonomy=False, include_deleted=False):
    """Prepare/check local files. Metadata downloads are managed by the user."""
    paths = [summary_path(key) for key in SUMMARY_URLS]
    if include_taxonomy:
        source = prepare_taxonomy(include_deleted=include_deleted)
        paths += [source / name for name in TAXONOMY_FILES + (("delnodes.dmp",) if include_deleted else ())]
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise RuntimeError("Missing local NCBI metadata:\n  " + "\n  ".join(map(str, missing)) +
                           "\nDownload missing metadata as described in 0.README.md, then rerun ./1.library.sh. "
                           "A local taxdump.tar.gz is extracted automatically when taxonomy files are missing. "
                           "The pipeline does not download metadata.")
    downloading = [path for path in paths if path.with_name(path.name + ".aria2").exists()]
    if downloading:
        raise RuntimeError("Metadata download is incomplete (.aria2 control file present):\n  " +
                           "\n  ".join(map(str, downloading)) + "\nFinish the aria2c download before building.")
    print(f"PASS: {len(paths)} local NCBI metadata files present", flush=True)

def parse_summary(data: bytes, accessions=None):
    if not data.endswith(b"\n"):
        raise ValueError("summary is empty or has an incomplete final line")
    header = None
    out = {}
    records = 0
    for line_number, line in enumerate(io.StringIO(data.decode("utf-8", errors="replace")), 1):
        raw = line.rstrip("\r\n")
        comment = raw.lstrip("#").strip()
        if raw.startswith("#") and comment.split("\t", 1)[0] == "assembly_accession":
            header = comment.split("\t")
            if not {"assembly_accession", "taxid"}.issubset(header):
                raise ValueError("summary header is missing assembly_accession/taxid")
            continue
        if not raw or raw.startswith("#"):
            continue
        vals = raw.split("\t")
        if not header or len(vals) != len(header):
            raise ValueError(f"invalid summary header or truncated row at line {line_number}")
        row = dict(zip(header, vals))
        if not re.fullmatch(r"GC[AF]_\d+\.\d+", row["assembly_accession"]) or not row["taxid"].isdigit():
            raise ValueError(f"invalid accession/taxid at line {line_number}")
        records += 1
        if accessions is None or row["assembly_accession"] in accessions:
            out[row["assembly_accession"]] = row
    if not header or not records:
        raise ValueError("summary has no assembly records")
    return out

def check_references(file_sizes=None, scan=None):
    assert_reference_policy()
    print(f"[references] Scanning {DIR_REF}", flush=True)
    sizes = {} if file_sizes is None else file_sizes
    unfinished = []

    def progress(group, count, state):
        message = "scanning..." if state == "start" else f"{count:,} FASTA files ({state})"
        print(f"[references] {group}: {message}", flush=True)

    if scan is None:
        files = reference_files(sizes=sizes, incomplete=unfinished, progress=progress)
    else:
        files, scanned_sizes, unfinished = scan
        sizes.update(scanned_sizes)
    if not files:
        raise SystemExit(f"No reference FASTA found under {DIR_REF}")

    missing_groups = [g for g in REF_GROUPS if not any(p.is_relative_to(DIR_REF/g) for p in files)]
    if missing_groups:
        raise SystemExit("Missing reference groups: " + ",".join(missing_groups))

    empty = [p for p in files if sizes[p] == 0]
    if empty:
        raise RuntimeError(f"Empty reference FASTA: {empty[0]} ({len(empty)} empty files)")
    if unfinished:
        raise RuntimeError(f"Reference download still incomplete: {unfinished[0]}; finish downloads first")
    print(f"PASS: {len(files):,} local reference files (presence/size check)", flush=True)
    return files


def load_metadata(files):
    metadata = {}
    receipts = {}
    loaded = set()
    accessions = {accession_from_filename(path) for path in files} - {None}
    for key, url in SUMMARY_URLS.items():
        source = summary_path(key)
        print(f"[metadata] Read local {source}", flush=True)
        raw = read_metadata_source(source)
        try:
            summary = parse_summary(raw, accessions=accessions)
        except ValueError as exc:
            raise RuntimeError(f"Cannot parse local metadata {source}: {exc}; download it again manually") from exc
        receipts[key] = {"url":url, "local_path":str(source), "sha256":hashlib.sha256(raw).hexdigest()}
        for row in summary.values():
            row["metadata_source"] = str(source)
            row["summary_group"] = key
        metadata.update(summary)
        print(f"[metadata] {key}: matched {len(summary):,} current assembly versions", flush=True)
        loaded.add(source)

    # Older selected mammal summaries kept the two human assemblies separately.
    human = summary_path("human")
    if human.is_file():
        print(f"[metadata] Read local human supplement {human}", flush=True)
        raw = read_metadata_source(human)
        for acc, row in parse_summary(raw, accessions=accessions).items():
            row["metadata_source"] = str(human)
            row["summary_group"] = "vertebrate_mammalian"
            metadata.setdefault(acc, row)
        receipts["human"] = {"local_path":str(human), "sha256":hashlib.sha256(raw).hexdigest()}
        loaded.add(human)

    # Selection snapshots may differ by a few assemblies. Fill only missing
    # accessions from other existing local snapshots, keeping primary records.
    wanted = {accession_from_filename(path) for path in files} - {None} - metadata.keys()
    if wanted:
        alternatives = {path for group in SUMMARY_URLS for path in summary_candidates(group)
                        if path not in loaded and path.is_file() and
                        not path.with_name(path.name + ".aria2").exists()}
        for source in sorted(alternatives, key=lambda path: (path.stat().st_size, str(path))):
            if not wanted:
                break
            print(f"[metadata] Check local supplement {source}", flush=True)
            raw = read_metadata_source(source)
            try:
                supplement = parse_summary(raw, accessions=wanted)
            except ValueError as exc:
                print(f"[metadata] Skip unusable optional supplement {source}: {exc}", flush=True)
                continue
            matched = wanted.intersection(supplement)
            if matched:
                for acc in matched:
                    supplement[acc]["metadata_source"] = str(source)
                    # Filename layouts differ across legacy summary snapshots.
                    supplement[acc]["summary_group"] = next(
                        group for group in SUMMARY_URLS if source in summary_candidates(group))
                    metadata[acc] = supplement[acc]
                wanted.difference_update(matched)
                receipts[f"supplement_{len(receipts)}"] = {
                    "local_path": str(source), "sha256": hashlib.sha256(raw).hexdigest(), "matched": len(matched)}
                print(f"[metadata] Matched {len(matched)} additional references locally", flush=True)

    return metadata, receipts


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def record_taxid_problems(rows, missing, resolved, source):
    problems = []
    groups = set()
    for row in rows:
        taxid = str(row.get("taxid", ""))
        species = str(row.get("species_taxid", ""))
        resolved_taxid = resolved.get(taxid, taxid)
        resolved_species = resolved.get(species, species)
        absent = {resolved_taxid, resolved_species} & missing
        if not absent:
            continue
        problems.append({**row, "resolved_taxid": resolved_taxid, "resolved_species_taxid": resolved_species,
                         "missing_taxids": ",".join(sorted(absent, key=int))})
        group = row.get("summary_group")
        if group in SUMMARY_URLS:
            groups.add(group)
    out = DIR_LIB/"reference"
    report = out/"taxid-unresolved.tsv"
    fields = ["assembly_accession", "filename", "organism_name", "taxid", "species_taxid",
              "resolved_taxid", "resolved_species_taxid", "missing_taxids", "metadata_source", "summary_group"]
    write_tsv(report, fields, problems)
    atomic_json(out/"taxid-validation.json", {
        "status": "needs_metadata", "taxonomy_source": str(source),
        "missing_taxids": sorted(missing, key=int), "references": problems,
    })
    # The pipeline generates download instructions, but never runs a downloader.
    urlfile = out/"metadata.url"
    lines = ["# Refresh assembly-to-taxid metadata; genome FASTA is retained.",
             "# Download with aria2c in WSL, then rerun ./1.library.sh."]
    for group in sorted(groups or SUMMARY_URLS.keys()):
        lines += [SUMMARY_URLS[group], f"  dir={DIR_META/'ncbi'}",
                  f"  out=assembly_summary.{group}.txt"]
    lines += ["https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
              f"  dir={DIR_META/'ncbi/taxonomy'}", "  out=taxdump.tar.gz"]
    urlfile.write_text("\n".join(lines) + "\n")
    return report, urlfile, len(problems)


def validate_taxonomy(rows, deleted_policy="error"):
    if deleted_policy not in ("error", "exclude"):
        raise ValueError("DELETED_TAXID_POLICY must be error or exclude")
    source = taxonomy_path()
    required = {str(row[key]) for row in rows for key in ("taxid", "species_taxid")
                if str(row.get(key, "")).isdigit() and int(row[key]) > 0}
    merged = {}
    with (source/"merged.dmp").open("rb", buffering=8 << 20) as handle:
        for line in handle:
            fields = line.split(b"\t|")
            if len(fields) >= 2:
                merged[fields[0].strip().decode("ascii")] = fields[1].strip().decode("ascii")
    resolved = {}
    for original in required:
        taxid, seen = original, set()
        while taxid in merged:
            if taxid in seen:
                raise RuntimeError(f"Cycle in local taxonomy mergers at taxid {taxid}")
            seen.add(taxid)
            taxid = merged[taxid]
        resolved[original] = taxid
    taxids = set(resolved.values())
    missing_by_file = {}
    for name in ("nodes.dmp", "names.dmp"):
        missing = set(taxids)
        print(f"[metadata] Check taxids in {source/name}", flush=True)
        # Large reads avoid thousands of small WSL/Windows filesystem requests.
        with (source/name).open("rb", buffering=8 << 20) as handle:
            for line in handle:
                fields = line.split(b"\t|")
                taxid = fields[0].strip().decode("ascii")
                if name == "names.dmp" and (len(fields) < 4 or fields[3].strip() != b"scientific name"):
                    continue
                missing.discard(taxid)
                if not missing:
                    break
        missing_by_file[name] = missing
    missing = set().union(*missing_by_file.values())
    exclusions = []
    if missing and deleted_policy == "exclude":
        deleted = set()
        with (source/"delnodes.dmp").open("rb", buffering=8 << 20) as handle:
            for line in handle:
                taxid = line.split(b"\t|", 1)[0].strip().decode("ascii")
                if taxid in missing:
                    deleted.add(taxid)
        kept = []
        for row in rows:
            bad = {resolved.get(str(row.get(key, "")), str(row.get(key, "")))
                   for key in ("taxid", "species_taxid")} & missing
            if bad and bad.issubset(deleted):
                exclusions.append({**row, "missing_taxids": ",".join(sorted(bad, key=int)),
                                   "reason": "deleted_taxid_without_valid_successor",
                                   "taxonomy_source": str(source)})
            else:
                kept.append(row)
        rows[:] = kept
        taxids = {resolved[str(row[key])] for row in rows for key in ("taxid", "species_taxid")
                  if str(row.get(key, "")) in resolved}
        missing &= taxids
    out = DIR_LIB/"reference"
    exclusion_report = out/"taxid-excluded.tsv"
    exclusion_fields = ["assembly_accession", "relative_path", "filename", "organism_name", "taxid", "species_taxid",
                        "missing_taxids", "reason", "metadata_source", "summary_group", "taxonomy_source"]
    write_tsv(exclusion_report, exclusion_fields, exclusions)
    if exclusions:
        print(f"[metadata] EXCLUDE {len(exclusions)} references with confirmed deleted taxids from database inputs; "
              f"source FASTA retained. List: {exclusion_report}", flush=True)
    if not rows:
        raise RuntimeError(f"No usable references remain after taxonomy validation; see {exclusion_report}")
    if missing:
        report, urlfile, count = record_taxid_problems(rows, missing, resolved, source)
        sample = ", ".join(sorted(missing, key=int)[:10])
        raise RuntimeError(f"Local taxonomy lacks {len(missing)} required taxids affecting {count} reference files "
                           f"(examples: {sample}).\nUnresolved mappings: {report}\n"
                           f"Refresh the assembly summaries with aria2c -c -i {urlfile}, then rerun ./1.library.sh.\n"
                           "Only an exact assembly accession/version with a valid taxid, or a merged.dmp successor, "
                           "is accepted. Missing IDs not confirmed in delnodes.dmp require further review; "
                           "no reference is assigned a guessed taxid.")

    # Keep an audit of official merge successors and successful metadata refreshes.
    prior = {}
    previous = out/"taxid-validation.json"
    if previous.is_file():
        try:
            saved = json.loads(previous.read_text())
            if isinstance(saved, dict):
                prior = {row["assembly_accession"]: row for row in saved.get("references", [])
                         if isinstance(row, dict) and row.get("assembly_accession")}
        except (ValueError, TypeError):
            pass  # An old diagnostic file is never used to determine a mapping.
    repairs = []
    for row in rows:
        original = dict(row)
        for key in ("taxid", "species_taxid"):
            if str(row.get(key, "")) in resolved:
                row[key] = resolved[str(row[key])]
        before = prior.get(row.get("assembly_accession"), original)
        if any(str(before.get(key, "")) != str(row.get(key, "")) for key in ("taxid", "species_taxid")):
            repairs.append({"assembly_accession": row.get("assembly_accession", ""),
                            "old_taxid": before.get("taxid", ""), "new_taxid": row.get("taxid", ""),
                            "old_species_taxid": before.get("species_taxid", ""),
                            "new_species_taxid": row.get("species_taxid", ""),
                            "metadata_source": row.get("metadata_source", ""), "taxonomy_source": str(source)})
    if repairs:
        write_tsv(out/"taxid-corrections.tsv", list(repairs[0]), repairs)
        print(f"[metadata] Repaired {len(repairs)} mappings; audit: {out/'taxid-corrections.tsv'}", flush=True)
    atomic_json(previous, {"status": "complete", "taxonomy_source": str(source),
                           "taxids_checked": len(taxids), "repairs": repairs,
                           "deleted_taxid_policy": deleted_policy, "active_references": len(rows),
                           "excluded_references": len(exclusions), "exclusions": exclusions})
    (out/"taxid-unresolved.tsv").unlink(missing_ok=True)
    print(f"PASS: {len(taxids):,} reference/species taxids exist in local taxonomy and scientific names", flush=True)
    return {"source": str(source), "taxids_checked": len(taxids),
            "merged_taxids": {old: new for old, new in resolved.items() if old != new},
            "repaired_references": len(repairs), "excluded_references": len(exclusions),
            "exclusion_report": str(exclusion_report), "deleted_taxid_policy": deleted_policy}


def main(include_taxonomy=False, deleted_policy="error", reference_scan=None):
    if deleted_policy not in ("error", "exclude"):
        raise ValueError("DELETED_TAXID_POLICY must be error or exclude")
    check_metadata(include_taxonomy=include_taxonomy, include_deleted=include_taxonomy and deleted_policy == "exclude")
    file_sizes = {}
    files = check_references(file_sizes=file_sizes, scan=reference_scan)
    metadata, receipts = load_metadata(files)
    outdir = DIR_LIB / "reference"
    outdir.mkdir(parents=True, exist_ok=True)

    rows, missing = [], []
    for path in files:
        acc = accession_from_filename(path)
        if not acc or acc not in metadata:
            missing.append((str(path), acc or ""))
            continue
        row = metadata[acc]
        rel = path.relative_to(DIR_REF)
        rows.append({
            "relative_path": str(rel),
            "filename": path.name,
            "assembly_accession": acc,
            "taxid": row.get("taxid",""),
            "species_taxid": row.get("species_taxid",""),
            "organism_name": row.get("organism_name",""),
            "assembly_level": row.get("assembly_level",""),
            "refseq_category": row.get("refseq_category",""),
            "bytes": file_sizes[path] if path in file_sizes else path.stat().st_size,
            "source_group": rel.parts[0],
            "metadata_source": row.get("metadata_source", ""),
            "summary_group": row.get("summary_group", ""),
        })

    if missing:
        with (outdir/"unmatched.tsv").open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f,delimiter="\t"); w.writerow(["path","assembly_accession"]); w.writerows(missing)
        raise SystemExit(f"{len(missing)} reference files not found in current RefSeq summaries; see {outdir/'unmatched.tsv'}")

    taxonomy = validate_taxonomy(rows, deleted_policy=deleted_policy) if include_taxonomy else None

    fields = ["relative_path", "filename", "assembly_accession", "taxid", "species_taxid",
              "organism_name", "assembly_level", "refseq_category", "bytes", "source_group", "metadata_source", "summary_group"]
    write_tsv(outdir/"inventory.tsv", fields, rows)

    # MetaCache assembly_summary.txt: first column MUST be filename, some column taxid.
    with (outdir/"assembly_summary.txt").open("w",encoding="utf-8") as f:
        f.write("# metag generated mapping from downloaded RefSeq assemblies\n")
        f.write("# assembly_accession\ttaxid\n")
        for r in rows:
            f.write(f"{r['filename']}\t{r['taxid']}\n")

    atomic_json(outdir/"reference.receipt.json", {
        "reference_root": str(DIR_REF),
        "files": len(rows),
        "source_files": len(files),
        "compressed_bytes": sum(int(r["bytes"]) for r in rows),
        "summary_sources": receipts,
        "taxonomy_validation": taxonomy,
        "policy": "reference FASTA retained permanently; no automatic cleanup",
    })
    print(f"PASS: {len(rows):,} active references from {len(files):,} source files; mapping={outdir/'assembly_summary.txt'}")
    return 0

if __name__=="__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
