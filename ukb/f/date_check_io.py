#!/usr/bin/env python3
"""Streaming local UKB date-audit input helpers. No network or third-party packages.

`linked` retains every matching raw GP/diagnosis record (including undated
records); `fields`/`fod` project wide TSV exports without loading the cohort.
Read/CTV3 codes are case-sensitive and their literal dots are significant.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import sys
import tempfile

csv.field_size_limit(1 << 30)
SOURCES = ("icd10", "icd9", "opcs4", "srd", "gp2", "gp3")
SPECIAL_DATES = {"1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"}
MISSING = {"", "NA", "N/A", "NaN", "nan", "NULL", "null"}
EVENT_FIELDS = ("eid", "source", "code", "date_raw", "origin", "granularity", "match_local", "match_official")


def open_text(path: Path):
    return (gzip.open if str(path).lower().endswith(".gz") else open)(path, "rt", encoding="utf-8-sig", newline="")


@contextlib.contextmanager
def atomic_text(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(fd)
    try:
        op = gzip.open if str(path).lower().endswith(".gz") else open
        with op(tmp, "wt", encoding="utf-8", newline="") as out:
            yield out
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def iter_rows(path: Path):
    """Read canonical UKB TSV; fail visibly if a row loses a field."""
    with open_text(path) as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header in {path}")
        reader.fieldnames = [x.strip() for x in reader.fieldnames]
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Duplicate column names in {path}")
        for rowno, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed TSV row {rowno} in {path}")
            yield rowno, row


def header(path: Path):
    with open_text(path) as inp:
        names = next(csv.reader(inp, delimiter="\t"), [])
    names = [x.strip() for x in names]
    if not names or len(set(names)) != len(names):
        raise ValueError(f"Missing/duplicate header in {path}")
    return names


def need_column(names, candidates, path, label):
    for name in candidates:
        if name in names:
            return name
    raise ValueError(f"{path}: missing {label} column (expected one of {', '.join(candidates)})")


def file_info(path: Path | None, source: str, status="available", digest=False):
    info = {"source": source, "path": str(path) if path else "", "status": status, "exists": bool(path and path.is_file())}
    if info["exists"]:
        st = path.stat()
        info.update(size_bytes=st.st_size, modified_utc=dt.datetime.fromtimestamp(st.st_mtime, dt.timezone.utc).isoformat())
        if digest:
            sha = hashlib.sha256()
            with path.open("rb") as inp:
                for block in iter(lambda: inp.read(1024 * 1024), b""):
                    sha.update(block)
            info["sha256"] = sha.hexdigest()
    return info


def choose_path(value, root: Path, stem: str):
    if str(value).lower() == "none":
        return None, "disabled"
    if str(value).lower() != "auto":
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f"Explicit input does not exist: {path}")
        return path, "available"
    for folder in (root / "rap" / "raw", root / "raw", root / "rap", root):
        for name in (stem + ".tsv.gz", stem + ".tsv", stem + ".tab.gz", stem + ".tab"):
            path = folder / name
            if path.is_file():
                return path.resolve(), "available"
    return None, "not_found"


def read_list(path: Path):
    if not path.is_file():
        return []
    out = []
    with open_text(path) as inp:
        for rowno, line in enumerate(inp, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) < 2:
                raise ValueError(f"{path}:{rowno}: expected two tab-separated columns")
            a, b = cells[0].strip(), cells[1].strip()
            if a and b:
                out.append((a, b))
    return out


def normal_code(source, code):
    code = str(code).strip()
    if source in ("gp2", "gp3"):
        return code
    if source == "srd":
        return re.sub(r"\.0+$", "", code)
    return code.upper().replace(".", "")


def field_id(name):
    match = re.fullmatch(r"p?(\d+)(?:_i\d+(?:_a\d+)?|_a\d+)?", str(name).strip())
    return int(match.group(1)) if match else None


class Definitions:
    def __init__(self, trait: str, common: Path, official_maps: Path):
        self.trait = trait
        self.inputs = []
        self.local_patterns = {}
        self.local_codes = {}
        self.matchers = {}
        for source in SOURCES:
            path = common / (source + ".lst")
            rows = read_list(path)
            pats = [a for a, b in rows if b == trait]
            self.local_patterns[source] = pats
            self.local_codes[source] = set()
            if source in ("gp2", "gp3"):
                self.local_codes[source] = {x.strip() for p in pats for x in p.split(",") if x.strip()}
                self.matchers[source] = []
            else:
                clean = [p.upper().replace(".", "") if source != "srd" else p for p in pats]
                # ICD/OPCS are prefix definitions. SRD is an exact integer-code
                # definition; 1471 must never accidentally match 14710.
                self.matchers[source] = [re.compile(r"^(?:" + p + r")" + (r"$" if source == "srd" else "")) for p in clean]
            self.inputs.append(file_info(path if path.is_file() else None, source + "_definitions", "available" if path.is_file() else "not_found", digest=True))
        fod_path, each_path = common / "fod.lst", common / "fod.each.lst"
        refs = [re.compile(r"^(?:" + a + r")$") for a, b in read_list(fod_path) if b == trait]
        each = read_list(each_path)
        self.ref_fields = sorted({"p" + str(field_id(a)) for a, _ in each if field_id(a) is not None and any(p.search(a) for p in refs)}, key=lambda x: int(x[1:]))
        self.ref_icd = sorted({normal_code("icd10", b)[:3] for a, b in each if "p" + str(field_id(a)) in self.ref_fields})
        if not self.ref_fields and not any(self.local_patterns.values()):
            raise ValueError(f"Trait {trait!r} is not defined in the available common/*.lst files")
        self.inputs.extend([file_info(p if p.is_file() else None, name, "available" if p.is_file() else "not_found", digest=True) for p, name in ((fod_path, "ref_definitions"), (each_path, "ref_field_dictionary"))])
        if not official_maps.is_file():
            raise ValueError(f"Official mapping file does not exist: {official_maps}")
        self.inputs.append(file_info(official_maps.resolve(), "official_maps", digest=True))
        self.official = {s: set() for s in SOURCES}
        self.dictionary = {s: set() for s in SOURCES}
        self.map_sources = set()
        with open_text(official_maps) as inp:
            reader = csv.DictReader(inp, delimiter="\t")
            if not reader.fieldnames or not {"source", "code", "icd10"} <= set(reader.fieldnames):
                raise ValueError("Official maps require source, code, icd10 columns")
            for row in reader:
                source = row["source"].strip().lower()
                if source not in self.dictionary:
                    continue
                code = normal_code(source, row["code"])
                if not code:
                    continue
                self.map_sources.add(source)
                self.dictionary[source].add(code)
                if normal_code("icd10", row["icd10"])[:3] in self.ref_icd:
                    self.official[source].add(code)
        self.official["icd10"] = set(self.ref_icd)
        self.dictionary["icd10"] = {normal_code("icd10", b)[:3] for _, b in each}
        self.map_sources.add("icd10")

    def local_defined(self, source):
        return bool(self.local_patterns.get(source))

    def official_defined(self, source):
        if source == "icd9" and not self.official[source]:
            # UKB coding 1836 is partial: it has no I48/4273 mapping.  Absence
            # of an official target mapping is unknown, never an official 0.
            return False
        return bool(self.ref_icd) and source in self.map_sources and source != "opcs4"

    def local_match(self, source, code):
        code = normal_code(source, code)
        if source in ("gp2", "gp3"):
            return code in self.local_codes[source]
        return any(p.search(code) is not None for p in self.matchers.get(source, []))

    def official_match(self, source, code):
        if not self.official_defined(source):
            return None
        code = normal_code(source, code)
        if source == "icd10":
            return code[:3] in self.official[source]
        # ICD9 mappings are explicit code-to-code rows, not an invented broad
        # ICD9 prefix mapping. Read2/CTV3/SRD also require exact dictionary codes.
        return code in self.official[source]

    def flags(self, source, code):
        local = self.local_match(source, code) if self.local_defined(source) else None
        official = self.official_match(source, code)
        return local, official

    def as_dict(self):
        return {"Y": self.trait, "ref_fields": self.ref_fields, "ref_icd10": self.ref_icd,
                "sources": {s: {"local_patterns": self.local_patterns[s], "local_codes": sorted(self.local_codes[s]),
                    "official_codes": sorted(self.official[s]), "local_defined": self.local_defined(s),
                    "official_defined": self.official_defined(s)} for s in SOURCES}}

    def comparison(self):
        out = []
        for source in SOURCES:
            local_defined, official_defined = self.local_defined(source), self.official_defined(source)
            if source in ("gp2", "gp3"):
                local = self.local_codes[source]
                scope = "Exact local Read/CTV3 code list versus bundled UKB mapping codes for the target three-character ICD10 group(s); case and literal dots preserved."
            else:
                local = {code for code in self.dictionary[source] if self.local_match(source, code)}
                scope = "Local pattern evaluated only on codes in the bundled official mapping dictionary; not an exhaustive clinical-code universe."
                if source == "icd10":
                    scope = "Comparison at UKB FOD three-character ICD10 level only; narrower local four-character definitions can differ within the same category and are not fully described by these sets."
            official = self.official[source]
            comparable = local_defined and official_defined
            out.append({"source": source, "local_defined": local_defined, "official_defined": official_defined,
                "local_only": sorted(local - official) if comparable else [],
                "official_only": sorted(official - local) if comparable else [],
                "shared": sorted(local & official) if comparable else [],
                "scope": scope if comparable else "One definition is unavailable or this source is not used in UKB first occurrences; absence is not a negative comparison."})
        return out

    def write(self, path: Path):
        with atomic_text(path) as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(("source", "kind", "value"))
            for source in SOURCES:
                for pattern in self.local_patterns[source]:
                    writer.writerow((source, "local_pattern", pattern))
                for code in sorted(self.local_codes[source]):
                    writer.writerow((source, "local_code", code))
                for code in sorted(self.official[source]):
                    writer.writerow((source, "official_code", code))
                writer.writerow((source, "local_defined", int(self.local_defined(source))))
                writer.writerow((source, "official_defined", int(self.official_defined(source))))
            for field in self.ref_fields:
                writer.writerow(("ref", "ref_field", field))
            for code in self.ref_icd:
                writer.writerow(("ref", "ref_icd", code))


def valid_calendar_date(value):
    value = str(value).strip()
    if value in MISSING:
        return None
    for fmt in ("%Y-%m-%d", "%d/%m/%Y", "%Y/%m/%d", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"):
        try:
            date = dt.datetime.strptime(value, fmt).date()
            return None if date.isoformat() in SPECIAL_DATES else date
        except ValueError:
            continue
    return None


def episode_date(row):
    fields = ("epistart", "admidate", "epiend", "disdate")
    for field in fields:
        if field in row and valid_calendar_date(row[field]) is not None:
            return row[field].strip(), field
    # Retain a special/invalid raw value for audit if no real date is available.
    for field in fields:
        if field in row and row[field].strip() not in MISSING:
            return row[field].strip(), field + "_invalid"
    return "", "missing_date"


def event_row(eid, source, code, raw_date, origin, granularity, local, official):
    return (eid, source, code, raw_date, origin, granularity, "" if local is None else int(local), "" if official is None else int(official))


def stream_gp(path: Path, defs: Definitions, writer, info):
    names = header(path)
    eid_col = need_column(names, ("eid", "participant.eid"), path, "eid")
    date_col = need_column(names, ("event_dt", "event_date", "date"), path, "event date")
    code_cols = {s: c for s, c in (("gp2", "read_2"), ("gp3", "read_3")) if c in names}
    if not code_cols:
        raise ValueError(f"{path}: neither read_2 nor read_3 is present")
    scanned = matched = output = 0
    for rowno, row in iter_rows(path):
        scanned += 1
        row_matched = False
        for source, col in code_cols.items():
            code = row[col].strip()
            if code in MISSING:
                continue
            local, official = defs.flags(source, code)
            if local or official:
                eid = row[eid_col].strip()
                if eid in MISSING:
                    raise ValueError(f"Matched GP record has no eid: {path}:{rowno}")
                origin = f"gp_clinical:row={rowno}"
                if "data_provider" in row:
                    origin += ";provider=" + row["data_provider"].strip()
                writer.writerow(event_row(eid, source, code, row[date_col].strip(), origin, "record", local, official))
                row_matched = True
                output += 1
        matched += int(row_matched)
    info.update(rows=scanned, matched_records=matched, output_rows=output, code_columns=list(code_cols.values()))


def stream_hes(diag: Path, hesin: Path | None, defs: Definitions, writer, infos, tmpdir: Path):
    names = header(diag)
    eid_col = need_column(names, ("eid", "participant.eid"), diag, "eid")
    index_col = need_column(names, ("ins_index",), diag, "episode identifier")
    code_cols = {s: c for s, c in (("icd10", "diag_icd10"), ("icd9", "diag_icd9")) if c in names}
    if not code_cols:
        raise ValueError(f"{diag}: neither diag_icd10 nor diag_icd9 is present")
    fd, dbpath = tempfile.mkstemp(prefix=".date-check-hes-", suffix=".sqlite", dir=tmpdir)
    os.close(fd)
    db = sqlite3.connect(dbpath)
    try:
        db.execute("PRAGMA journal_mode=OFF")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("CREATE TABLE diagnoses (eid TEXT, ins_index TEXT, source TEXT, code TEXT, ml INTEGER, mo INTEGER, origin TEXT, seen INTEGER DEFAULT 0)")
        batch = []
        scanned = matched_rows = stored = 0
        for rowno, row in iter_rows(diag):
            scanned += 1
            row_matched = False
            for source, col in code_cols.items():
                code = row[col].strip()
                if code in MISSING:
                    continue
                local, official = defs.flags(source, code)
                if local or official:
                    eid = row[eid_col].strip()
                    if eid in MISSING:
                        raise ValueError(f"Matched diagnosis has no eid: {diag}:{rowno}")
                    origin = f"hesin_diag:row={rowno};ins_index=" + row[index_col].strip()
                    if "arr_index" in row:
                        origin += ";arr_index=" + row["arr_index"].strip()
                    batch.append((eid, row[index_col].strip(), source, code, None if local is None else int(local), None if official is None else int(official), origin))
                    row_matched = True
                    stored += 1
            matched_rows += int(row_matched)
            if len(batch) >= 5000:
                db.executemany("INSERT INTO diagnoses(eid,ins_index,source,code,ml,mo,origin) VALUES (?,?,?,?,?,?,?)", batch)
                batch.clear()
        if batch:
            db.executemany("INSERT INTO diagnoses(eid,ins_index,source,code,ml,mo,origin) VALUES (?,?,?,?,?,?,?)", batch)
        db.execute("CREATE INDEX episode_key ON diagnoses(eid,ins_index)")
        db.commit()
        infos["hesin_diag"].update(rows=scanned, matched_records=matched_rows, matched_codes=stored, code_columns=list(code_cols.values()))
        output = matched_episodes = scanned_episodes = 0
        if hesin:
            hnames = header(hesin)
            heid = need_column(hnames, ("eid", "participant.eid"), hesin, "eid")
            hindex = need_column(hnames, ("ins_index",), hesin, "episode identifier")
            if not any(x in hnames for x in ("epistart", "admidate", "epiend", "disdate")):
                raise ValueError(f"{hesin}: no supported episode date columns")
            for rowno, row in iter_rows(hesin):
                scanned_episodes += 1
                key = (row[heid].strip(), row[hindex].strip())
                # Missing identifiers must not accidentally join two unknowns.
                if key[0] in MISSING or key[1] in MISSING:
                    continue
                hits = db.execute("SELECT eid,source,code,ml,mo,origin FROM diagnoses WHERE eid=? AND ins_index=?", key).fetchall()
                if not hits:
                    continue
                matched_episodes += 1
                date, chosen = episode_date(row)
                for eid, source, code, ml, mo, origin in hits:
                    origin += f";hesin_row={rowno};date_field={chosen}"
                    writer.writerow(event_row(eid, source, code, date, origin, "episode", ml, mo))
                    output += 1
                db.execute("UPDATE diagnoses SET seen=1 WHERE eid=? AND ins_index=?", key)
            infos["hesin"].update(rows=scanned_episodes, matched_records=matched_episodes)
        missing = 0
        for eid, source, code, ml, mo, origin in db.execute("SELECT eid,source,code,ml,mo,origin FROM diagnoses WHERE seen=0"):
            writer.writerow(event_row(eid, source, code, "", origin + ";episode_record_unavailable", "episode", ml, mo))
            missing += 1
            output += 1
        infos["hesin_diag"].update(output_rows=output, unmatched_episode_codes=missing)
    finally:
        db.close()
        for name in (dbpath, dbpath + "-journal", dbpath + "-wal", dbpath + "-shm"):
            if os.path.exists(name):
                os.unlink(name)


def linked(args):
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    root = Path(args.phenotype_dir).expanduser().resolve()
    common = Path(args.common_dir).expanduser().resolve()
    defs = Definitions(args.Y, common, Path(args.official_maps).expanduser().resolve())
    defs.write(outdir / "definitions.tsv")
    gp, gp_status = choose_path(args.gp_records, root, "gp_clinical")
    hesin, hesin_status = choose_path(args.hesin, root, "hesin")
    diag, diag_status = choose_path(args.hesin_diag, root, "hesin_diag")
    infos = {"gp_records": file_info(gp, "gp_records", gp_status), "hesin": file_info(hesin, "hesin", hesin_status), "hesin_diag": file_info(diag, "hesin_diag", diag_status)}
    with atomic_text(outdir / "linked-events.tsv") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(EVENT_FIELDS)
        if gp:
            stream_gp(gp, defs, writer, infos["gp_records"])
        if diag:
            stream_hes(diag, hesin, defs, writer, infos, outdir)
    limitations = [
        "A missing input or undefined code set is unavailable evidence, not evidence that participants are disease-free.",
        "Raw GP dates and unmatched/undated hospital diagnoses are retained; only downstream date validation can determine usable event dates.",
        "One GP record may contain both Read2 and CTV3 codes; both code matches are retained, so counts are not independent episode counts.",
        "Hospital episode date priority: first valid epistart, admidate, epiend, disdate. The last two are fallback approximations and are identified in origin; an all-invalid raw date is retained for audit.",
        "The official mapping file is a versioned reconstruction aid. Local data releases, source coverage, date imputation and UKB mapping versions may differ from the published first-occurrence field.",
    ]
    if not defs.ref_icd:
        limitations.append("No exact UKB FOD field mapping for this trait was found in fod.lst plus fod.each.lst; official reconstruction cannot be evaluated.")
    if defs.ref_icd and not defs.official["icd9"]:
        limitations.append("The published partial ICD9 mapping contains no codes for target " + ", ".join(defs.ref_icd) + "; official ICD9 reconstruction is unavailable. Local codes (including 4273 for I48) are never silently substituted into the official candidate.")
    if not gp:
        limitations.append("Raw GP records unavailable/disabled: only any separately loaded summary/AOD data can describe GP dates.")
    if not diag:
        limitations.append("Raw HES diagnoses unavailable/disabled: hospital summary AOD contains first dates per distinct code, not a complete admission history.")
    elif not hesin:
        limitations.append("HES diagnoses available but episodes unavailable/disabled: matching diagnosis rows are retained without an episode date.")
    metadata = {
        "Y": args.Y, "definitions": defs.as_dict(), "mapping_comparison": defs.comparison(),
        "sources": list(infos.values()), "inputs": defs.inputs + list(infos.values()), "limitations": limitations,
        "raw_gp2_available": "read_2" in infos["gp_records"].get("code_columns", []),
        "raw_gp3_available": "read_3" in infos["gp_records"].get("code_columns", []),
        "raw_icd10_available": "diag_icd10" in infos["hesin_diag"].get("code_columns", []),
        "raw_icd9_available": "diag_icd9" in infos["hesin_diag"].get("code_columns", []),
        "raw_hospital_episode_dates_available": bool(diag and hesin),
        "first_occurrence_maps_version": "UKB Resource 593 mappings (January 2022 GP/SRD release) plus the partial UKB coding 1836 ICD9 mapping retrieved 2026-10-30; exact bundled file SHA256 is recorded in inputs.",
    }
    with atomic_text(outdir / "metadata.tsv") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(("key", "value"))
        for key, value in metadata.items():
            writer.writerow((key, json.dumps(value, ensure_ascii=False, separators=(",", ":"))))
    print(json.dumps({"Y": args.Y, "output_dir": str(outdir), "sources": list(infos.values())}, ensure_ascii=False))


def project_fields(input_path: Path, selected, outputs):
    """Write field projections in one pass, always including eid.

    A set selects every instance of its field IDs. A dictionary maps field IDs
    to an all-instances flag; false selects only baseline/singular columns.
    """
    paths = [p.resolve() for p in outputs]
    if input_path.resolve() in paths or len(set(paths)) != len(paths):
        raise ValueError("Input and output paths must be distinct")
    names = header(input_path)
    eid = need_column(names, ("eid", "participant.eid"), input_path, "eid")
    def selected_column(name, ids):
        fid = field_id(name)
        if fid not in ids:
            return False
        if not isinstance(ids, dict) or ids[fid]:
            return True
        return re.fullmatch(r"p?\d+(?:_i0(?:_a0)?|_a0)?", name) is not None
    columns = [[eid] + [n for n in names if n != eid and selected_column(n, ids)] for ids in selected]
    if any(len(cols) == 1 for cols in columns) and len(outputs) == 1:
        raise ValueError("None of the requested fields/instances is present")
    count = 0
    with contextlib.ExitStack() as stack:
        writers = [csv.writer(stack.enter_context(atomic_text(p)), delimiter="\t", lineterminator="\n") for p in outputs]
        for writer, cols in zip(writers, columns):
            writer.writerow(["eid" if c == eid else c for c in cols])
        indices = [[names.index(c) for c in cols] for cols in columns]
        # Avoid a dictionary with tens of thousands of phenotype keys for
        # every participant; only retain the requested column indices.
        with open_text(input_path) as inp:
            reader = csv.reader(inp, delimiter="\t")
            next(reader)
            for rowno, row in enumerate(reader, 2):
                if not row:
                    continue
                if len(row) != len(names):
                    raise ValueError(f"Malformed TSV row {rowno} in {input_path}")
                for writer, ids in zip(writers, indices):
                    writer.writerow([row[i] for i in ids])
                count += 1
    return {"input": str(input_path), "rows": count, "outputs": [{"path": str(p), "columns": cols} for p, cols in zip(outputs, columns)]}


def fields(args):
    path = Path(args.input).expanduser().resolve()
    mapping = getattr(args, "mapping", None)
    if mapping:
        # Dataset lists use field ID, label, and an optional third-column 'all'.
        # Keep the established baseline policy for every other third-column value.
        ids = {}
        mapping_path = Path(mapping).expanduser().resolve()
        with open_text(mapping_path) as inp:
            for rowno, line in enumerate(inp, 1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                cells = line.rstrip("\r\n").split("\t")
                if cells[0].strip() == "eid":
                    continue
                fid = field_id(cells[0])
                if fid is None:
                    raise ValueError(f"{mapping_path}:{rowno}: invalid UKB field ID {cells[0]!r}")
                all_instances = len(cells) >= 3 and cells[2].strip().lower() == "all"
                ids[fid] = ids.get(fid, False) or all_instances
    else:
        ids = {int(x.strip().removeprefix("p")) for x in args.fields.split(",") if x.strip()}
    if not ids or not any(field_id(n) in ids for n in header(path)):
        raise ValueError("None of the requested fields is present")
    print(json.dumps(project_fields(path, [ids], [Path(args.output)]), ensure_ascii=False))


def fod(args):
    path = Path(args.input).expanduser().resolve()
    ids = {field_id(a) for a, _ in read_list(Path(args.mapping)) if field_id(a) is not None}
    if not ids or any(i % 2 for i in ids):
        raise ValueError("FOD mapping must contain even-numbered date field identifiers")
    if not any(field_id(n) in ids for n in header(path)):
        raise ValueError("No selected FOD date fields are present in the export")
    print(json.dumps(project_fields(path, [ids, {i + 1 for i in ids}], [Path(args.date_output), Path(args.source_output)]), ensure_ascii=False))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("linked", help="Stream matching GP/HES records and audit code definitions")
    for flag in ("Y", "phenotype-dir", "common-dir", "output-dir", "official-maps"):
        p.add_argument("--" + flag, required=True)
    for flag in ("gp-records", "hesin", "hesin-diag"):
        p.add_argument("--" + flag, default="auto")
    p.set_defaults(func=linked)
    p = sub.add_parser("fields", help="Extract exact base fields from a wide phenotype TSV")
    for flag in ("input", "output"):
        p.add_argument("--" + flag, required=True)
    select = p.add_mutually_exclusive_group(required=True)
    select.add_argument("--fields", help="Comma-separated base field IDs; include all instances")
    select.add_argument("--mapping", help="Dataset list: field ID, label, optional 'all'; otherwise baseline/singular only")
    p.set_defaults(func=fields)
    p = sub.add_parser("fod", help="Extract FOD dates and paired provenance fields in one pass")
    for flag in ("input", "mapping", "date-output", "source-output"):
        p.add_argument("--" + flag, required=True)
    p.set_defaults(func=fod)
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except (OSError, ValueError, csv.Error, re.error, sqlite3.Error) as exc:
        parser.exit(2, f"date_check_io: {exc}\n")


if __name__ == "__main__":
    main()
