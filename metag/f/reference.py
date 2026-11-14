"""Frozen reference preparation and an uncapped native Kraken2 HPC database.

Input FASTA files are never rewritten or removed. Tasks retain complete records;
their size target refers to compressed input bytes, not a chromosome split size.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import tarfile
import threading
import urllib.request
from contextlib import contextmanager

from common import atomic_json, check_space, file_state, fingerprint, locked, sha256, utcnow

SCHEMA = "metag-hpc-reference-v1"
DB_SCHEMA = "metag-hpc-database-v1"
SOURCE_COMMIT = "8c190b1b668825935dbf6dee5f969227dc8269bb"
GIB = 1024 ** 3
GROUPS = ("viral", "archaea", "bacteria", "fungi", "protozoa", "plant",
          "invertebrate", "vertebrate_mammalian", "vertebrate_other")
TAX_FILES = ("nodes.dmp", "names.dmp", "merged.dmp")
ACC = re.compile(r"^(GC[AF]_\d+\.\d+)(?:_|\.|$)")
FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fas")
DNA = frozenset(b"ACGTURYSWKMBDHVNXacgturyswkmbdhvnx")
LOWER_TO_X = bytes.maketrans(b"abcdefghijklmnopqrstuvwxyz", b"x" * 26)
MAX_MASK_RECORD = (1 << 31) - 1  # Pinned k2mask stores single-record positions in int32_t.
PARAMETERS = {"k": 35, "l": 31, "seed": "1" * 17 + "01" * 7,
              "load_factor": 0.7, "block_size": 16384, "subblock_size_rule": "ceil(block_size / threads)",
              "min_taxid_bits": 0, "maximum_capacity": None,
              "minimum_acceptable_hash_value": 0, "deterministic": True}


def _paths(cfg):
    lib = Path(cfg["library_dir"])
    meta = Path(cfg["metadata_dir"])
    return {"lib": lib, "meta": meta, "raw": Path(cfg["reference_raw"]).resolve(),
            "tax": meta / "taxonomy", "db": lib / "db",
            "prepared": lib / "prepared", "plan": lib / "reference_plan.json",
            "lock": (Path(cfg["work_dir"]) / "jobs" / "reference.lock") if cfg.get("storage_layout") == "separated" else lib / "reference.lock"}


def _read_json(path):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Invalid or missing checkpoint {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected an object in checkpoint {path}")
    return value


def _sync_file(handle):
    handle.flush()
    os.fsync(handle.fileno())


def _sync_dir(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _owned_stage(path, owner):
    """Recreate only this workflow's interrupted private stage, never a foreign dir."""
    path = Path(path)
    if path.is_symlink():
        raise RuntimeError(f"Refusing a symlink staging directory: {path}")
    if path.exists():
        marker = path / "owner.json"
        if not path.is_dir() or not marker.is_file() or marker.is_symlink() or _read_json(marker) != owner:
            raise RuntimeError(f"Refusing an unowned/incompatible staging directory: {path}")
        shutil.rmtree(path)
    path.mkdir()
    atomic_json(path / "owner.json", owner)


def _remove_owned_stage(path, owner):
    if path.is_symlink() or not path.is_dir() or _read_json(path / "owner.json") != owner:
        raise RuntimeError(f"Staging ownership changed before cleanup: {path}")
    shutil.rmtree(path)


def _file_record(path, *, checksum=True):
    before = file_state(path)
    digest = sha256(path) if checksum else None
    if file_state(path) != before:
        raise RuntimeError(f"File changed while hashing: {path}")
    return {"state": before, **({"sha256": digest} if checksum else {})}


def _check_record(path, record, *, full=False):
    if not Path(path).is_file() or file_state(path) != record["state"]:
        raise RuntimeError(f"Frozen file changed or disappeared: {path}")
    if full and "sha256" in record and sha256(path) != record["sha256"]:
        raise RuntimeError(f"Checksum mismatch: {path}")


def _software(cfg):
    folder = Path(cfg["tools_dir"]) / "kraken"
    receipt = _read_json(folder / "install.json")
    if receipt.get("source_tag") != "2.17.2" or receipt.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError("Kraken installation is not the pinned official 2.17.2 commit; run setup")
    names = ("k2mask", "estimate_capacity", "build_db", "classify")
    files = {}
    for name in names:
        path = folder / "bin" / name
        if not path.is_file() or not os.access(path, os.X_OK):
            raise RuntimeError(f"Missing native executable {path}; run setup")
        digest = sha256(path)
        if receipt.get("files", {}).get(name) != digest:
            raise RuntimeError(f"Native executable does not match its installation receipt: {path}")
        files[name] = digest
    return {"source_tag": "2.17.2", "source_commit": SOURCE_COMMIT, "files": files}


def _inventory(cfg):
    root = _paths(cfg)["raw"]
    if not root.is_dir():
        raise RuntimeError(f"Reference input directory does not exist: {root}")
    result = []
    for folder, dirs, names in os.walk(root, followlinks=False):
        for name in dirs:
            if name.startswith(".rsync-partial"):
                raise RuntimeError(f"Rsync staging directory inside reference input: {Path(folder)/name}; move it outside reference_raw after the upload finishes")
            if (Path(folder) / name).is_symlink():
                raise RuntimeError(f"Reference directory symlink requires an explicit real input path: {Path(folder)/name}")
        for name in names:
            path = Path(folder) / name
            if name.endswith((".aria2", ".partial", ".part", ".filepart", ".tmp")):
                raise RuntimeError(f"Unfinished reference download/file: {path}")
            stem = name[:-3] if name.endswith(".gz") else name
            if not stem.lower().endswith(FASTA_SUFFIXES) and stem not in ("UniVec", "UniVec_Core"):
                if name.lower().endswith((".md5", ".sha256", ".md5sum", ".txt", ".tsv", ".url")):
                    continue
                raise RuntimeError(f"Unsupported file in reference input: {path}; expected DNA FASTA[.gz]")
            if path.is_symlink() or not path.is_file():
                raise RuntimeError(f"Reference must be a regular, non-symlink file: {path}")
            state = file_state(path)
            if state["size"] == 0:
                raise RuntimeError(f"Empty reference file: {path}")
            match = ACC.match(name)
            result.append({"path": path.relative_to(root).as_posix(), "state": state,
                           "accession": match.group(1) if match else None})
    result.sort(key=lambda x: x["path"])
    if not result:
        raise RuntimeError(f"No DNA FASTA files found below {root}")
    return result


def _custom_map(cfg, inventory):
    value = cfg.get("reference_taxid_map")
    if not value:
        return {}, None
    path = Path(value)
    root = _paths(cfg)["raw"].resolve()
    allowed = {row["path"] for row in inventory}
    result = {}
    record = _file_record(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"path", "taxid"}.issubset(reader.fieldnames):
            raise RuntimeError(f"Custom taxid TSV needs path and taxid columns: {path}")
        for row in reader:
            candidate = Path(row["path"])
            candidate = candidate if candidate.is_absolute() else root / candidate
            try:
                relative = candidate.resolve().relative_to(root).as_posix()
                taxid = int(row["taxid"])
            except (ValueError, TypeError) as exc:
                raise RuntimeError(f"Invalid custom mapping row: {row}") from exc
            if relative not in allowed or taxid <= 0 or relative in result:
                raise RuntimeError(f"Unknown/duplicate path or invalid taxid in custom mapping: {row}")
            result[relative] = taxid
    _check_record(path, record)
    return result, {"path": str(path), **record}


def _taxonomy_ids(folder):
    parents, merged = {}, {}
    with (folder / "nodes.dmp").open() as handle:
        for number, line in enumerate(handle, 1):
            cols = [x.strip() for x in line.split("|")]
            try:
                node, parent = int(cols[0]), int(cols[1])
            except (ValueError, IndexError) as exc:
                raise RuntimeError(f"Invalid taxonomy node at {folder}/nodes.dmp:{number}") from exc
            if node in parents or node <= 0 or parent <= 0:
                raise RuntimeError(f"Duplicate/invalid taxonomy node {node}")
            parents[node] = parent
    if parents.get(1) != 1 or any(parent not in parents for parent in parents.values()):
        raise RuntimeError("Taxonomy must have root 1 and a parent for every node")
    with (folder / "merged.dmp").open() as handle:
        for line in handle:
            if line.strip():
                cols = line.split("|")
                merged[int(cols[0])] = int(cols[1])
    if not (folder / "names.dmp").stat().st_size:
        raise RuntimeError("Taxonomy names.dmp is empty")

    def resolve(taxid):
        seen = set()
        while taxid in merged:
            if taxid in seen:
                raise RuntimeError(f"Cycle in merged taxonomy for {taxid}")
            seen.add(taxid)
            taxid = merged[taxid]
        if taxid not in parents:
            raise RuntimeError(f"Taxid {taxid} is absent/deleted in the frozen NCBI taxonomy; update metadata or provide an explicit mapping")
        node, seen = taxid, set()
        while node != 1:
            if node in seen:
                raise RuntimeError(f"Taxonomy parent cycle for {taxid}")
            seen.add(node)
            node = parents[node]
        return taxid
    return resolve


def _summary_rows(path, wanted):
    result, header, records = {}, None, 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.endswith("\n"):
                raise RuntimeError(f"Incomplete assembly summary line at {path}:{number}")
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                fields = line.lstrip("#").strip().split("\t")
                if fields[0] == "assembly_accession":
                    header = fields
                continue
            if not line:
                continue
            fields = line.split("\t")
            if not header or "taxid" not in header or len(fields) != len(header):
                raise RuntimeError(f"Malformed assembly summary at {path}:{number}")
            accession = fields[0]
            taxid = fields[header.index("taxid")]
            if not re.fullmatch(r"GC[AF]_\d+\.\d+", accession) or not taxid.isdigit():
                raise RuntimeError(f"Invalid accession/taxid at {path}:{number}")
            records += 1
            if accession in wanted:
                value = int(taxid)
                if accession in result and result[accession] != value:
                    raise RuntimeError(f"Conflicting taxids for {accession} in {path}")
                result[accession] = value
    if not header or not records:
        raise RuntimeError(f"Assembly summary has no records: {path}")
    return result


def _download(url, destination):
    print(f"DOWNLOAD metadata {url}", flush=True)
    request = urllib.request.Request(url, headers={"User-Agent": "metag-hpc/1 (reference metadata)"})
    temporary = destination.with_name(destination.name + ".download")
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as out:
            shutil.copyfileobj(response, out, 4 * 1024 ** 2)
            _sync_file(out)
        if not temporary.stat().st_size:
            raise RuntimeError(f"Empty metadata response: {url}")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _tax_bundle(source, destination):
    if not source or not Path(source).exists():
        return False
    source = Path(source)
    candidates = [source, source / "taxonomy", source / "ncbi/taxonomy"]
    for folder in candidates:
        if all((folder / name).is_file() for name in TAX_FILES):
            records = {name: file_state(folder / name) for name in TAX_FILES + ("delnodes.dmp",)
                       if (folder / name).is_file()}
            for name in records:
                if (folder / (name + ".aria2")).exists():
                    raise RuntimeError(f"Incomplete taxonomy download: {folder/name}")
            for name in TAX_FILES + ("delnodes.dmp",):
                if (folder / name).is_file():
                    shutil.copyfile(folder / name, destination / name)
            if any(file_state(folder / name) != state for name, state in records.items()):
                raise RuntimeError(f"Taxonomy changed during import: {folder}")
            return True
    archives = [source / "taxdump.tar.gz", source / "taxonomy/taxdump.tar.gz",
                source / "ncbi/taxonomy/taxdump.tar.gz", source / "ncbi/taxdump.tar.gz"]
    for archive in archives:
        if archive.is_file():
            if archive.with_name(archive.name + ".aria2").exists():
                raise RuntimeError(f"Incomplete taxonomy archive download: {archive}")
            before = file_state(archive)
            _extract_taxdump(archive, destination)
            if file_state(archive) != before:
                raise RuntimeError(f"Taxonomy archive changed during extraction: {archive}")
            return True
    return False


def _extract_taxdump(archive, destination):
    with tarfile.open(archive, "r:gz") as bundle:
        members = {member.name.removeprefix("./"): member for member in bundle.getmembers()}
        for name in TAX_FILES + ("delnodes.dmp",):
            member = members.get(name)
            if member is None:
                if name in TAX_FILES:
                    raise RuntimeError(f"Taxdump archive lacks {name}: {archive}")
                continue
            if not member.isfile():
                raise RuntimeError(f"Taxdump member is not a regular file: {name}")
            with bundle.extractfile(member) as inp, (destination / name).open("wb") as out:
                shutil.copyfileobj(inp, out)
                _sync_file(out)


def _supplement_taxonomy(folder, xml_path):
    """Add only missing taxa from a saved NCBI efetch response; never guess ancestry."""
    import xml.etree.ElementTree as ET
    xml_path = Path(xml_path)
    if not xml_path.is_file():
        return
    nodes_path = folder / "nodes.dmp"
    records = {}
    for line in nodes_path.read_text().splitlines():
        cols = [v.strip() for v in line.split("|")]
        records[int(cols[0])] = cols
    taxa = ET.parse(xml_path).getroot().findall("Taxon")
    added = {}
    for taxon in taxa:
        chain = taxon.findall("LineageEx/Taxon") + [taxon]
        parent = 1
        for item in chain:
            taxid = int(item.findtext("TaxId"))
            actual_parent = int(item.findtext("ParentTaxId") or parent)
            if taxid in records and int(records[taxid][1]) != actual_parent:
                raise RuntimeError(f"NCBI supplement conflicts with taxonomy parent for {taxid}")
            if taxid not in records:
                if item is not taxon or actual_parent not in records:
                    raise RuntimeError(f"NCBI supplement needs an unavailable ancestor for {taxid}")
                name, rank = item.findtext("ScientificName"), item.findtext("Rank")
                if not name or not rank or any(c in name for c in "\t\n|"):
                    raise RuntimeError(f"Malformed NCBI supplement for {taxid}")
                # Inherit division/genetic-code fields from the verified parent.
                fields = list(records[actual_parent][:13])
                fields[0:3] = [str(taxid), str(actual_parent), rank]
                fields[3] = ""
                fields[12] = ""
                if taxon.findtext("GeneticCode/GCId"):
                    fields[6] = taxon.findtext("GeneticCode/GCId")
                    fields[7] = "0"
                if taxon.findtext("MitoGeneticCode/MGCId"):
                    fields[8] = taxon.findtext("MitoGeneticCode/MGCId")
                    fields[9] = "0"
                records[taxid] = fields
                added[taxid] = (fields, name)
            parent = taxid
    with nodes_path.open("a") as out, (folder / "names.dmp").open("a") as names:
        for taxid, (fields, name) in sorted(added.items()):
            out.write("\t|\t".join(fields) + "\t|\n")
            names.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    deleted = folder / "delnodes.dmp"
    if deleted.exists() and added:
        lines = deleted.read_text().splitlines(True)
        deleted.write_text("".join(line for line in lines if int(line.split("|")[0]) not in added))
    shutil.copyfile(xml_path, folder / "ncbi_supplement.xml")
    print(f"Imported {len(added)} missing taxonomy nodes from saved NCBI efetch", flush=True)


def _metadata_snapshot(cfg):
    paths = _paths(cfg)
    receipt = _read_json(paths["meta"] / "metadata.json")
    if receipt.get("schema") != SCHEMA or receipt.get("status") != "complete" or \
            fingerprint(receipt.get("context")) != receipt.get("identity"):
        raise RuntimeError("Metadata checkpoint is incomplete or has an invalid identity")
    for name, record in receipt["context"]["files"].items():
        _check_record(paths["meta"] / name, record)
    custom = receipt["context"].get("custom_map")
    if (str(cfg.get("reference_taxid_map")) if cfg.get("reference_taxid_map") else None) != \
            (custom["path"] if custom else None):
        raise RuntimeError("Custom reference mapping differs from the frozen metadata")
    if custom:
        _check_record(Path(custom["path"]), custom)
    return receipt


def metadata(cfg, source=None, download=True):
    """Import/download only metadata; resolve every raw file without guessing taxa."""
    paths = _paths(cfg)
    with locked(paths["lock"]):
        inventory = _inventory(cfg)
        custom, custom_record = _custom_map(cfg, inventory)
        inventory_id = fingerprint(inventory)
        marker = paths["meta"] / "metadata.json"
        if marker.is_file():
            try:
                existing = _metadata_snapshot(cfg)
                if existing["context"]["raw_inventory"] == inventory_id and \
                        existing["context"]["custom_map"] == custom_record:
                    print("REUSE complete reference metadata", flush=True)
                    return existing
            except RuntimeError:
                if paths["plan"].exists():
                    raise
        if paths["plan"].exists() or (paths["db"] / "build-state.json").exists():
            raise RuntimeError("Reference inputs/metadata are frozen by an existing plan; use a separate library for changed inputs")
        paths["meta"].mkdir(parents=True, exist_ok=True)
        check_space(paths["meta"], float(cfg.get("min_free_gib", 100)))
        if source is None:
            candidates = [Path(cfg["root"]) / "annot/metag/0.meta", paths["raw"].parent / "0.meta"]
            source = next((path for path in candidates if path.is_dir()), None)
        source = Path(source) if source is not None else None
        if source is not None and not source.is_dir():
            raise RuntimeError(f"Metadata import directory does not exist: {source}")
        stage = paths["meta"] / ".bootstrap.partial"
        owner = {"schema": SCHEMA, "kind": "metadata-bootstrap", "reference_raw": str(paths["raw"]),
                 "metadata_dir": str(paths["meta"]), "source": str(source) if source else None}
        _owned_stage(stage, owner)
        (stage / "taxonomy").mkdir()
        (stage / "summaries").mkdir()
        if not (_tax_bundle(paths["meta"], stage / "taxonomy") or _tax_bundle(source, stage / "taxonomy")):
            if not download:
                raise RuntimeError("Missing nodes.dmp/names.dmp/merged.dmp or taxdump.tar.gz; import metadata or enable official NCBI download")
            archive = stage / "taxdump.tar.gz"
            _download("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", archive)
            _extract_taxdump(archive, stage / "taxonomy")
            archive.unlink()
        _supplement_taxonomy(stage / "taxonomy", paths["meta"] / "taxonomy_supplement.xml")
        resolve = _taxonomy_ids(stage / "taxonomy")
        # Match bundled Kraken scripts/download_genomic_library.sh: UniVec is other sequences.
        for row in inventory:
            name = Path(row["path"]).name.removesuffix(".gz")
            if name in ("UniVec", "UniVec_Core"):
                custom.setdefault(row["path"], 28384)
        for path in custom:
            custom[path] = resolve(custom[path])
        missing_nonassembly = [row["path"] for row in inventory if not row["accession"] and row["path"] not in custom]
        if missing_nonassembly:
            raise RuntimeError("Non-assembly FASTA requires explicit reference_taxid_map path/taxid TSV; unmapped: " + ", ".join(missing_nonassembly[:10]))
        wanted = {row["accession"] for row in inventory if row["path"] not in custom}
        wanted.discard(None)
        matches, used = {}, []
        imported = set()
        # The explicitly requested snapshot wins; existing managed summaries fill gaps.
        for origin in (source, paths["meta"] / "summaries"):
            if origin and origin.is_dir():
                for candidate in sorted(origin.rglob("*assembly_summary*")):
                    if stage in candidate.parents:
                        continue
                    if not candidate.is_file() or candidate.name.endswith((".aria2", ".partial", ".download")):
                        continue
                    name = candidate.name
                    if name in ("assembly_summary.txt", "assembly_summary.txt.gz") and candidate.parent.name in GROUPS:
                        name = f"assembly_summary.{candidate.parent.name}.txt" + (".gz" if name.endswith(".gz") else "")
                    if not name.endswith((".txt", ".txt.gz")) or name in imported:
                        continue
                    target = stage / "summaries" / name
                    if candidate.with_name(candidate.name + ".aria2").exists():
                        raise RuntimeError(f"Incomplete metadata download: {candidate}")
                    before = file_state(candidate)
                    shutil.copyfile(candidate, target)
                    if file_state(candidate) != before:
                        raise RuntimeError(f"Assembly summary changed during import: {candidate}")
                    imported.add(name)
        def consume(path):
            found = _summary_rows(path, wanted - matches.keys())
            for accession, taxid in found.items():
                matches[accession] = (resolve(taxid), path.name)
            used.append(path.name)
        summaries = sorted((stage / "summaries").iterdir(),
                           key=lambda p: ("historical" in p.name, "genbank" in p.name, p.name))
        for path in summaries:
            if wanted - matches.keys():
                consume(path)
        if download:
            for group in GROUPS:
                if not wanted - matches.keys():
                    break
                name = f"assembly_summary.{group}.txt"
                path = stage / "summaries" / name
                if name not in used:
                    if not path.exists():
                        _download(f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{group}/assembly_summary.txt", path)
                    consume(path)
            for historical in (False, True):
                for namespace, prefix in (("refseq", "GCF_"), ("genbank", "GCA_")):
                    if not any(acc.startswith(prefix) for acc in wanted - matches.keys()):
                        continue
                    suffix = "_historical" if historical else ""
                    name = f"assembly_summary_{namespace}{suffix}.txt"
                    path = stage / "summaries" / name
                    if name not in used:
                        if not path.exists():
                            _download("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/" + name, path)
                        consume(path)
        missing = sorted(wanted - matches.keys())
        if missing:
            raise RuntimeError(f"{len(missing)} assembly versions have no taxid mapping; no raw reference is excluded. Import current/historical RefSeq/GenBank summaries or provide reference_taxid_map. Examples: {', '.join(missing[:15])}")
        mapping = stage / "reference_taxids.tsv"
        with mapping.open("w", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(("path", "accession", "taxid", "source"))
            for row in inventory:
                taxid, origin = (custom[row["path"]], "custom") if row["path"] in custom else matches[row["accession"]]
                writer.writerow((row["path"], row["accession"] or "", taxid, origin))
            _sync_file(out)
        if fingerprint(_inventory(cfg)) != inventory_id:
            raise RuntimeError("Reference inventory changed during metadata preparation")
        if custom_record:
            _check_record(Path(custom_record["path"]), custom_record)
        # There is no frozen plan at this point. Publish complete files, then the receipt.
        published = []
        for category in ("taxonomy", "summaries"):
            destination = paths["meta"] / category
            if destination.is_symlink():
                raise RuntimeError(f"Refusing a symlink managed metadata directory: {destination}")
            destination.mkdir(exist_ok=True)
            for candidate in (stage / category).iterdir():
                target = destination / candidate.name
                candidate.replace(target)
                with target.open("rb") as handle:
                    os.fsync(handle.fileno())
                published.append(target)
            _sync_dir(destination)
        mapping.replace(paths["meta"] / mapping.name)
        published.append(paths["meta"] / mapping.name)
        context = {"raw_inventory": inventory_id, "custom_map": custom_record,
                   "files": {str(path.relative_to(paths["meta"])): _file_record(path, checksum=False) for path in published}}
        result = {"schema": SCHEMA, "status": "complete", "identity": fingerprint(context),
                  "context": context, "references": len(inventory), "assemblies": len(wanted),
                  "custom_references": len(custom), "created_at": utcnow()}
        atomic_json(marker, result)
        _remove_owned_stage(stage, owner)
        return result


def _load_plan(cfg, verify_inventory=False):
    paths = _paths(cfg)
    plan = _read_json(paths["plan"])
    if plan.get("schema") != SCHEMA or plan.get("status") != "complete" or \
            fingerprint(plan.get("context")) != plan.get("identity"):
        raise RuntimeError("Invalid or incomplete reference plan")
    context = plan["context"]
    if plan.get("tasks") != context["tasks"]:
        raise RuntimeError("Reference plan task list differs from its fingerprint")
    if context["reference_raw"] != str(paths["raw"]) or \
            context["metadata_identity"] != _metadata_snapshot(cfg)["identity"]:
        raise RuntimeError("Reference plan does not match current raw directory/metadata")
    if context["masking"] != bool(cfg.get("masking", True)) or context["software"] != _software(cfg):
        raise RuntimeError("Reference masking/toolchain differs from the frozen plan")
    if verify_inventory and fingerprint(_inventory(cfg)) != context["raw_inventory"]:
        raise RuntimeError("Raw references changed after planning; use a separate library")
    return plan


def plan_references(cfg):
    paths = _paths(cfg)
    with locked(paths["lock"]):
        if paths["plan"].exists():
            return _load_plan(cfg, verify_inventory=True)
        meta = _metadata_snapshot(cfg)
        inventory = _inventory(cfg)
        if fingerprint(inventory) != meta["context"]["raw_inventory"]:
            raise RuntimeError("Raw inventory changed since metadata resolution; rerun metadata before planning")
        mapping = {}
        with (paths["meta"] / "reference_taxids.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                mapping[row["path"]] = int(row["taxid"])
        if set(mapping) != {row["path"] for row in inventory}:
            raise RuntimeError("Reference mapping must cover every and only current raw FASTA file")
        for index, row in enumerate(inventory):
            row.update(index=index, taxid=mapping[row["path"]])
        target = float(cfg.get("prepare_target_gib", 2)) * GIB
        max_files = int(cfg.get("prepare_max_files", 1000))
        if target <= 0 or max_files <= 0:
            raise RuntimeError("prepare_target_gib and prepare_max_files must be positive")
        tasks, indices, size = [], [], 0
        for row in inventory:
            if indices and (size + row["state"]["size"] > target or len(indices) >= max_files):
                tasks.append({"task_id": len(tasks) + 1, "files": indices, "compressed_bytes": size})
                indices, size = [], 0
            indices.append(row["index"])
            size += row["state"]["size"]
        if indices:
            tasks.append({"task_id": len(tasks) + 1, "files": indices, "compressed_bytes": size})
        # Array tasks read only their small manifest, not all reference stat records.
        manifests = paths["lib"] / "reference_manifests"
        manifests.mkdir(exist_ok=True)
        for task in tasks:
            manifest = manifests / f"task{task['task_id']:06d}.json"
            atomic_json(manifest, {"task_id": task["task_id"],
                                   "files": [inventory[index] for index in task["files"]]})
            task["manifest"] = {"path": str(manifest.relative_to(paths["lib"])), **_file_record(manifest)}
        context = {"reference_raw": str(paths["raw"]), "metadata_identity": meta["identity"],
                   "raw_inventory": meta["context"]["raw_inventory"], "masking": bool(cfg.get("masking", True)),
                   "software": _software(cfg), "prepare_target_bytes": int(target),
                   "prepare_max_files": max_files, "tasks": tasks,
                   "record_policy": "whole-record-unique-id-lowercase-to-x-v1"}
        plan = {"schema": SCHEMA, "status": "complete", "identity": fingerprint(context),
                "context": context, "tasks": tasks, "references": len(inventory),
                "compressed_bytes": sum(row["state"]["size"] for row in inventory), "created_at": utcnow()}
        atomic_json(paths["plan"], plan)
        return plan


def _task_rows(cfg, task):
    record = task["manifest"]
    path = _paths(cfg)["lib"] / record["path"]
    _check_record(path, record)
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise RuntimeError(f"Reference task manifest checksum mismatch: {path}")
    document = json.loads(data)
    if document.get("task_id") != task["task_id"] or \
            [row["index"] for row in document["files"]] != task["files"]:
        raise RuntimeError(f"Reference task manifest does not match plan: {path}")
    _check_record(path, record)
    return document["files"]


@contextmanager
def _interruptible():
    old = signal.getsignal(signal.SIGTERM)
    def stop(signum, frame):
        raise KeyboardInterrupt("Terminated")
    signal.signal(signal.SIGTERM, stop)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, old)


def _stop(process):
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()


def _emit_fasta(cfg, rows, destination, mapping, stats):
    headers = hashlib.sha256()
    total_records, total_bases, longest = 0, 0, 0
    sources = []
    for row in rows:
        path = _paths(cfg)["raw"] / row["path"]
        if file_state(path) != row["state"]:
            raise RuntimeError(f"Raw reference changed after planning: {path}")
        records, bases, record_bases = 0, 0, 0
        with path.open("rb", buffering=0) as raw:
            compressed = io.BufferedReader(raw, buffer_size=1024 ** 2)
            decoded = gzip.GzipFile(fileobj=compressed) if path.name.endswith(".gz") else compressed
            try:
                for line in decoded:
                    line = line.rstrip(b"\r\n")
                    if not line:
                        continue
                    if line.startswith(b">"):
                        if len(line) == 1 or (records and record_bases == 0):
                            raise RuntimeError(f"Empty FASTA header/record in {path}")
                        longest = max(longest, record_bases)
                        records += 1
                        record_bases = 0
                        seqid = f"metag_f{row['index']:010d}_r{records:010d}|kraken:taxid|{row['taxid']}"
                        header = b">" + seqid.encode() + b"\n"
                        destination.write(header)
                        headers.update(header)
                        mapping.write((seqid + "\t" + str(row["taxid"]) + "\n").encode())
                    else:
                        if not records or not set(line).issubset(DNA):
                            raise RuntimeError(f"Invalid DNA FASTA sequence in {path}, record {records}")
                        if bool(cfg.get("masking", True)) and record_bases + len(line) > MAX_MASK_RECORD:
                            raise RuntimeError(f"FASTA record exceeds the pinned k2mask signed-position limit: {path}, record {records}. No completed cache will be published; handle this record with a separately validated masking strategy.")
                        destination.write(line + b"\n")
                        bases += len(line)
                        record_bases += len(line)
                if not records or not record_bases:
                    raise RuntimeError(f"Empty/truncated FASTA in {path}")
                longest = max(longest, record_bases)
            finally:
                if decoded is not compressed:
                    decoded.close()
                compressed.close()
        if file_state(path) != row["state"]:
            raise RuntimeError(f"Raw reference changed during preparation: {path}")
        sources.append({"path": row["path"], "state": row["state"],
                        "records": records, "bases": bases})
        total_records += records
        total_bases += bases
    stats.update(records=total_records, bases=total_bases, longest_record=longest,
                 header_sha256=headers.hexdigest(), sources=sources)


def _prepare_native(cfg, rows, stage):
    stats, errors = {}, []
    output = stage / "sequences.fna"
    map_path = stage / "seqid2taxid.map"
    if not bool(cfg.get("masking", True)):
        with output.open("wb") as out, map_path.open("wb") as mapping:
            _emit_fasta(cfg, rows, out, mapping, stats)
            _sync_file(out)
            _sync_file(mapping)
        return stats
    threads = int(cfg.get("prepare_threads", 8))
    if threads <= 0:
        raise RuntimeError("prepare_threads must be positive")
    command = [str(Path(cfg["tools_dir"]) / "kraken/bin/k2mask"), "-threads",
               str(threads), "-outfmt", "fasta"]
    with (stage / "mask.log").open("w") as log, map_path.open("wb") as mapping, \
            output.open("wb") as out, _interruptible():
        log.write("COMMAND " + json.dumps(command) + "\n")
        log.flush()
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=log, start_new_session=True)
        def produce():
            try:
                _emit_fasta(cfg, rows, process.stdin, mapping, stats)
                process.stdin.close()
            except BaseException as exc:
                errors.append(exc)
                _stop(process)
        producer = threading.Thread(target=produce, daemon=True)
        producer.start()
        headers, bases, records = hashlib.sha256(), 0, 0
        try:
            for line in process.stdout:
                if line.startswith(b">"):
                    headers.update(line)
                    records += 1
                else:
                    line = line.translate(LOWER_TO_X)
                    bases += len(line.rstrip(b"\r\n"))
                out.write(line)
            returncode = process.wait()
            producer.join()
            if errors:
                raise RuntimeError(f"Reference input producer failed: {errors[0]}") from errors[0]
            if returncode:
                raise RuntimeError(f"k2mask failed with exit {returncode}; see {stage/'mask.log'}")
            if (records, bases, headers.hexdigest()) != (stats["records"], stats["bases"], stats["header_sha256"]):
                raise RuntimeError("k2mask changed reference headers/order or sequence lengths")
            _sync_file(mapping)
            _sync_file(out)
        finally:
            _stop(process)
            if process.stdin:
                try:
                    process.stdin.close()
                except OSError:
                    pass
            process.stdout.close()
            producer.join(timeout=10)
        log.write(f"EXIT_CODE {returncode}\n")
    return stats


def _prepared_receipt(cfg, plan, task):
    folder = _paths(cfg)["prepared"] / f"task{task['task_id']:06d}"
    receipt = _read_json(folder / "receipt.json")
    if receipt.get("schema") != SCHEMA or receipt.get("status") != "complete" or \
            receipt.get("plan_identity") != plan["identity"] or receipt.get("task_id") != task["task_id"] or \
            fingerprint(receipt.get("context")) != receipt.get("identity"):
        raise RuntimeError(f"Invalid or incompatible prepared task checkpoint: {folder}")
    if receipt["context"]["task"] != task:
        raise RuntimeError(f"Task definition differs from receipt: {folder}")
    for name, record in receipt["context"]["files"].items():
        _check_record(folder / name, record)
    for row in _task_rows(cfg, task):
        if file_state(_paths(cfg)["raw"] / row["path"]) != row["state"]:
            raise RuntimeError(f"Raw input changed: {row['path']}")
    return folder, receipt


def prepare_task(cfg, task_id):
    paths = _paths(cfg)
    with locked(paths["lock"], shared=True):
        plan = _load_plan(cfg)
        tasks = plan["context"]["tasks"]
        if not 1 <= int(task_id) <= len(tasks):
            raise RuntimeError(f"Prepare task ID must be in 1..{len(tasks)}")
        task = tasks[int(task_id) - 1]
        paths["prepared"].mkdir(parents=True, exist_ok=True)
        final = paths["prepared"] / f"task{int(task_id):06d}"
        stage = paths["prepared"] / f".task{int(task_id):06d}.partial"
        with locked(paths["prepared"] / f"task{int(task_id):06d}.lock"):
            if final.exists():
                _, receipt = _prepared_receipt(cfg, plan, task)
                print(f"REUSE prepare task {task_id}", flush=True)
                return receipt
            check_space(paths["prepared"], float(cfg.get("min_free_gib", 100)))
            owner = {"schema": SCHEMA, "kind": "prepare-task", "plan_identity": plan["identity"],
                     "task_id": task["task_id"]}
            _owned_stage(stage, owner)
            rows = _task_rows(cfg, task)
            stats = _prepare_native(cfg, rows, stage)
            context = {"plan_identity": plan["identity"], "task": task, "stats": stats,
                       "files": {name: _file_record(stage / name, checksum=False) for name in ("sequences.fna", "seqid2taxid.map")}}
            receipt = {"schema": SCHEMA, "status": "complete", "plan_identity": plan["identity"],
                       "task_id": task["task_id"], "identity": fingerprint(context), "context": context,
                       "execution": {"threads": int(cfg.get("prepare_threads", 8)), "masking": bool(cfg.get("masking", True))},
                       "created_at": utcnow()}
            atomic_json(stage / "receipt.json", receipt)
            _sync_dir(stage)
            stage.rename(final)
            _sync_dir(paths["prepared"])
            return receipt


def _run_feed(command, sources, log_path, stdout_path=None):
    """Feed complete prepared files while checking their recorded paths and sizes."""
    errors = []
    with log_path.open("w") as log, _interruptible():
        log.write("COMMAND " + json.dumps(command) + "\n")
        log.flush()
        output = stdout_path.open("wb") if stdout_path else log
        environment = dict(os.environ)
        if "-p" in command:
            environment["OMP_NUM_THREADS"] = command[command.index("-p") + 1]
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=output, stderr=log,
                                   start_new_session=True, env=environment)
        def produce():
            try:
                for path, record in sources:
                    _check_record(path, record)
                    digest = hashlib.sha256() if "sha256" in record else None
                    with path.open("rb") as inp:
                        for block in iter(lambda: inp.read(4 * 1024 ** 2), b""):
                            if digest is not None:
                                digest.update(block)
                            process.stdin.write(block)
                    _check_record(path, record)
                    if digest is not None and digest.hexdigest() != record["sha256"]:
                        raise RuntimeError(f"Prepared reference checksum mismatch: {path}")
                process.stdin.close()
            except BaseException as exc:
                errors.append(exc)
                _stop(process)
        producer = threading.Thread(target=produce, daemon=True)
        producer.start()
        try:
            returncode = process.wait()
            producer.join()
            if errors:
                raise RuntimeError(f"Prepared input stream failed: {errors[0]}") from errors[0]
            if returncode:
                raise RuntimeError(f"Native program failed with exit {returncode}: {command[0]}; see {log_path}")
        finally:
            _stop(process)
            try:
                process.stdin.close()
            except OSError:
                pass
            producer.join(timeout=10)
            if stdout_path:
                output.close()
        log.write(f"EXIT_CODE {returncode}\n")


def _native_opts(path):
    data = path.read_bytes()
    # Pinned native ABI on supported 64-bit little-endian Linux: size_t/u64,
    # bool + alignment, then the minimum accepted hash at byte offset 40.
    if len(data) < 52 or struct.calcsize("P") != 8:
        raise RuntimeError("Unsupported/truncated native Kraken2 options layout")
    k, l, _, _ = struct.unpack_from("<QQQQ", data)
    minimum = struct.unpack_from("<Q", data, 40)[0]
    if k != PARAMETERS["k"] or l != PARAMETERS["l"] or not data[32] or minimum != 0:
        raise RuntimeError("Built database options do not prove the requested uncapped DNA configuration")
    return {"k": k, "l": l, "minimum_acceptable_hash_value": minimum}


def _native_index_structure(folder, *, expected_capacity=None):
    """Check pinned native file layouts with two 32-byte reads, not a full scan.

    The bundled native hash writer does not check ostream write/close errors.
    A zero exit status and a nonempty hash file therefore do not establish that
    all cells reached disk. These checks also apply to legacy valid receipts;
    they neither rewrite those receipts nor change the build identity.
    """
    folder = Path(folder)
    headers, sizes = {}, {}
    for name in ("hash.k2d", "taxo.k2d"):
        path = folder / name
        if path.is_symlink() or not path.is_file():
            raise RuntimeError(f"Native index must be a regular non-symlink file: {path}")
        with path.open("rb") as handle:
            headers[name] = handle.read(32)
            sizes[name] = os.fstat(handle.fileno()).st_size
        if len(headers[name]) != 32:
            raise RuntimeError(f"Truncated native index header: {path}")

    capacity, used, key_bits, value_bits = struct.unpack("<QQQQ", headers["hash.k2d"])
    if not capacity or used > capacity or not 0 < value_bits < 32 or key_bits + value_bits != 32:
        raise RuntimeError(f"Invalid native hash header: {folder / 'hash.k2d'}")
    if expected_capacity is not None and (
            isinstance(expected_capacity, bool) or not isinstance(expected_capacity, int)
            or expected_capacity <= 0 or capacity != expected_capacity):
        raise RuntimeError(f"Native hash capacity {capacity} differs from the build estimate {expected_capacity}")
    expected_hash_bytes = 32 + capacity * 4
    if sizes["hash.k2d"] != expected_hash_bytes:
        raise RuntimeError(f"Native hash length mismatch: {folder / 'hash.k2d'}; "
                           f"expected {expected_hash_bytes} bytes, found {sizes['hash.k2d']}")

    # Pinned 64-bit little-endian ABI: 'K2TAXDAT', three size_t fields,
    # node_count records of seven uint64_t fields, then two string arrays.
    header = headers["taxo.k2d"]
    if header[:8] != b"K2TAXDAT":
        raise RuntimeError(f"Invalid native taxonomy magic: {folder / 'taxo.k2d'}")
    nodes, names_bytes, ranks_bytes = struct.unpack("<QQQ", header[8:])
    if nodes < 2 or not names_bytes or not ranks_bytes:
        raise RuntimeError(f"Invalid native taxonomy header: {folder / 'taxo.k2d'}")
    expected_taxonomy_bytes = 32 + nodes * 56 + names_bytes + ranks_bytes
    if sizes["taxo.k2d"] != expected_taxonomy_bytes:
        raise RuntimeError(f"Native taxonomy length mismatch: {folder / 'taxo.k2d'}; "
                           f"expected {expected_taxonomy_bytes} bytes, found {sizes['taxo.k2d']}")
    if nodes > (1 << value_bits):
        raise RuntimeError("Native hash taxid bits cannot encode the native taxonomy node count")
    return {"capacity": capacity, "used": used, "key_bits": key_bits,
            "value_bits": value_bits, "taxonomy_nodes": nodes,
            "hash_bytes": expected_hash_bytes, "taxonomy_bytes": expected_taxonomy_bytes}


def database_snapshot(cfg):
    """Fast read-only validation; never rehash the TB-scale database per sample."""
    paths = _paths(cfg)
    marker = _read_json(paths["db"] / "database.json")
    if marker.get("schema") != DB_SCHEMA or marker.get("status") != "complete" or \
            fingerprint(marker.get("context")) != marker.get("identity"):
        raise RuntimeError("Global database has no valid complete checkpoint")
    if marker.get("db_dir") != str(paths["db"]):
        raise RuntimeError("Database checkpoint belongs to a different directory")
    if set(marker.get("files", {})) != {"hash.k2d", "taxo.k2d", "opts.k2d"} or \
            not set(TAX_FILES).issubset(marker.get("taxonomy_files", {})):
        raise RuntimeError("Database checkpoint lacks required native/taxonomy files")
    if marker.get("taxonomy_dir") != str(paths["db"] / "taxonomy") or \
            marker.get("parameters") != marker["context"].get("parameters"):
        raise RuntimeError("Database taxonomy path/parameters differ from its checkpoint")
    for name, record in marker["files"].items():
        _check_record(paths["db"] / name, record)
    for name, record in marker["taxonomy_files"].items():
        _check_record(paths["db"] / "taxonomy" / name, record)
    if marker["files"] != marker["context"]["files"] or marker["taxonomy_files"] != marker["context"]["taxonomy_files"]:
        raise RuntimeError("Database output records do not match its fingerprint")
    _native_opts(paths["db"] / "opts.k2d")
    resource = marker["context"].get("resource_estimate", {})
    _native_index_structure(paths["db"], expected_capacity=resource.get("capacity"))
    return marker


def build_database(cfg):
    paths = _paths(cfg)
    with locked(paths["lock"], shared=True):
        plan = _load_plan(cfg, verify_inventory=True)
        paths["db"].mkdir(parents=True, exist_ok=True)
        with locked(paths["db"] / "build.lock"):
            if (paths["db"] / "database.json").exists():
                complete = database_snapshot(cfg)
                if complete["context"]["plan_identity"] != plan["identity"]:
                    raise RuntimeError("Completed database belongs to another reference plan")
                return complete
            state_path = paths["db"] / "build-state.json"
            allowed = {"build.lock", "build-state.json"}
            if not state_path.exists() and any(p.name not in allowed for p in paths["db"].iterdir()):
                raise RuntimeError("Refusing an unowned/incomplete database directory without build-state.json")
            prepared = [_prepared_receipt(cfg, plan, task) for task in plan["context"]["tasks"]]
            build_context = {"plan_identity": plan["identity"], "parameters": PARAMETERS,
                             "software": _software(cfg),
                             "tasks": [receipt["identity"] for _, receipt in prepared]}
            build_id = fingerprint(build_context)
            if state_path.exists() and _read_json(state_path).get("identity") != build_id:
                raise RuntimeError("Existing build state belongs to different inputs; use a separate library")
            atomic_json(state_path, {"schema": DB_SCHEMA, "status": "building", "identity": build_id,
                                     "context": build_context, "started_at": utcnow()})
            check_space(paths["db"], float(cfg.get("min_free_gib", 100)))
            taxonomy = paths["db"] / "taxonomy"
            taxonomy.mkdir(exist_ok=True)
            metadata_state = _metadata_snapshot(cfg)
            for name in TAX_FILES + ("delnodes.dmp",):
                source = paths["tax"] / name
                if source.exists():
                    record = metadata_state["context"]["files"].get("taxonomy/" + name)
                    if record is None:
                        continue
                    _check_record(source, record)
                    destination = taxonomy / name
                    if destination.exists():
                        if file_state(destination) != record["state"]:
                            raise RuntimeError(f"Frozen database taxonomy differs from metadata: {destination}")
                    else:
                        temporary = destination.with_name(name + ".partial")
                        shutil.copyfile(source, temporary)
                        with temporary.open("rb") as handle:
                            os.fsync(handle.fileno())
                        temporary.replace(destination)
            _sync_dir(taxonomy)
            frozen_taxonomy = {name: _file_record(taxonomy / name) for name in TAX_FILES + ("delnodes.dmp",)
                               if (taxonomy / name).is_file()}
            library = paths["db"] / "library"
            library.mkdir(exist_ok=True)
            sources = []
            mapping = paths["db"] / "seqid2taxid.map"
            temporary = mapping.with_name(mapping.name + ".partial")
            with temporary.open("wb") as out:
                for folder, receipt in prepared:
                    link = library / folder.name
                    if link.is_symlink():
                        if link.resolve() != folder.resolve():
                            raise RuntimeError(f"Foreign prepared-library link: {link}")
                    elif link.exists():
                        raise RuntimeError(f"Expected a link to prepared references: {link}")
                    else:
                        link.symlink_to(folder, target_is_directory=True)
                    record = receipt["context"]["files"]["seqid2taxid.map"]
                    _check_record(folder / "seqid2taxid.map", record, full=True)
                    with (folder / "seqid2taxid.map").open("rb") as inp:
                        shutil.copyfileobj(inp, out, 4 * 1024 ** 2)
                    sources.append((folder / "sequences.fna", receipt["context"]["files"]["sequences.fna"]))
                _sync_file(out)
            temporary.replace(mapping)
            frozen_mapping = _file_record(mapping)
            binary = Path(cfg["tools_dir"]) / "kraken/bin"
            threads = int(cfg.get("build_threads", 32))
            if threads <= 0:
                raise RuntimeError("build_threads must be positive")
            base = ["-k", "35", "-l", "31", "-S", PARAMETERS["seed"], "-p", str(threads)]
            capacity_path = paths["db"] / "capacity.json"
            if capacity_path.exists():
                estimate_receipt = _read_json(capacity_path)
                if estimate_receipt.get("build_identity") != build_id or estimate_receipt.get("status") != "complete" or \
                        estimate_receipt.get("identity") != fingerprint({"build_identity": build_id,
                            "distinct_minimizers_estimate": estimate_receipt.get("distinct_minimizers_estimate")}):
                    raise RuntimeError("Capacity estimate belongs to another/incomplete build")
                estimate = int(estimate_receipt["distinct_minimizers_estimate"])
            else:
                _run_feed([str(binary / "estimate_capacity"), *base], sources,
                          paths["db"] / "estimate.log", paths["db"] / "estimate.stdout")
                text = (paths["db"] / "estimate.stdout").read_text().strip()
                if not text.isdigit() or int(text) <= 0:
                    raise RuntimeError(f"Invalid native capacity estimate: {text!r}")
                estimate = int(text)
                estimate_context = {"build_identity": build_id, "distinct_minimizers_estimate": estimate}
                atomic_json(capacity_path, {"status": "complete", **estimate_context,
                                           "identity": fingerprint(estimate_context), "created_at": utcnow()})
            capacity = (estimate + 8192) * 10 // 7
            hash_bytes = capacity * 4
            longest = max(receipt["context"]["stats"]["longest_record"] for _, receipt in prepared)
            # Hash table plus a conservative working allowance, not a peak-memory guarantee.
            estimated_ram = math.ceil(hash_bytes * 1.15) + 8 * GIB + longest * (threads + 2)
            budget = float(cfg.get("build_memory_gib", 4096)) * GIB
            resource = {"capacity": capacity, "estimated_hash_bytes": hash_bytes,
                        "estimated_ram_bytes": estimated_ram, "budget_bytes": int(budget),
                        "no_capacity_cap": True}
            atomic_json(paths["db"] / "resource_estimate.json", resource)
            if estimated_ram > budget:
                raise RuntimeError(f"Uncapped build estimated RAM {estimated_ram/GIB:.2f} GiB exceeds budget {budget/GIB:.2f} GiB; request a larger-memory job. No downsampling was enabled.")
            check_space(paths["db"], float(cfg.get("min_free_gib", 100)) + hash_bytes * 1.01 / GIB)
            stage = paths["db"] / ".building.partial"
            owner = {"schema": DB_SCHEMA, "kind": "native-build", "build_identity": build_id}
            _owned_stage(stage, owner)
            command = [str(binary / "build_db"), *base,
                       "-H", str(stage / "hash.k2d"), "-t", str(stage / "taxo.k2d"),
                       "-o", str(stage / "opts.k2d"), "-n", str(taxonomy), "-m", str(mapping),
                       "-c", str(capacity), "-B", "16384", "-b", str(math.ceil(16384 / threads)), "-r", "0"]
            _run_feed(command, sources, paths["db"] / "build.log")
            native_options = _native_opts(stage / "opts.k2d")
            _native_index_structure(stage, expected_capacity=capacity)
            for name in ("hash.k2d", "taxo.k2d", "opts.k2d"):
                file = stage / name
                if not file.is_file() or file.stat().st_size == 0:
                    raise RuntimeError(f"Native build did not produce {name}")
                with file.open("rb") as handle:
                    os.fsync(handle.fileno())
                file.replace(paths["db"] / name)
            files = {name: _file_record(paths["db"] / name) for name in ("hash.k2d", "taxo.k2d", "opts.k2d")}
            for name, record in frozen_taxonomy.items():
                _check_record(taxonomy / name, record, full=True)
            _check_record(mapping, frozen_mapping, full=True)
            taxonomy_files = frozen_taxonomy
            # Confirm frozen task/raw/metadata states after the long build, before publishing.
            _load_plan(cfg, verify_inventory=True)
            for task in plan["context"]["tasks"]:
                _prepared_receipt(cfg, plan, task)
            context = {**build_context, "files": files, "taxonomy_files": taxonomy_files,
                       "resource_estimate": resource, "native_options": native_options}
            complete = {"schema": DB_SCHEMA, "status": "complete", "identity": fingerprint(context),
                        "context": context, "db_dir": str(paths["db"]), "files": files,
                        "taxonomy_dir": str(taxonomy), "taxonomy_files": taxonomy_files,
                        "parameters": dict(PARAMETERS), "completed_at": utcnow()}
            atomic_json(paths["db"] / "database.json", complete)
            _remove_owned_stage(stage, owner)
            _sync_dir(paths["db"])
            return complete
