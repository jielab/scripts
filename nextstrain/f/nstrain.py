#!/usr/bin/env python3
"""
Utility script for a small hChikV / GISAID / Nextstrain build.

The shell script should call this script for all data-curation tasks:
  1) prepare FASTA + metadata from GISAID exports
  2) add country centroids and write lat_longs.tsv
  3) write auspice_config.json
  4) count metadata states before optional augur traits steps

No country coordinates are hard-coded here. Coordinates are read from
countries.centroids.tsv, which should have four tab-separated columns:
  resolution  value  latitude  longitude
Example:
  country     China  35.86166  104.195397
"""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ACC_RE = re.compile(r"EPI_ISL_\d+")
UNKNOWN_VALUES = {"", "unknown", "unk", "na", "n/a", "nan", "none", "null", "not_applicable"}

# Small alias table only for common naming mismatches between GISAID and country centroid files.
# Coordinates themselves are still read from countries.centroids.tsv.
COUNTRY_ALIASES = {
    "USA": "United_States",
    "U.S.A.": "United_States",
    "US": "United_States",
    "United_States_of_America": "United_States",
    "United States of America": "United_States",
    "United States": "United_States",
    "UK": "United_Kingdom",
    "U.K.": "United_Kingdom",
    "United Kingdom": "United_Kingdom",
    "UAE": "United_Arab_Emirates",
    "United Arab Emirates": "United_Arab_Emirates",
    "South Korea": "South_Korea",
    "Republic of Korea": "South_Korea",
    "Korea, South": "South_Korea",
    "North Korea": "North_Korea",
    "Russian Federation": "Russia",
    "Czechia": "Czech_Republic",
    "Viet Nam": "Vietnam",
    "Côte d’Ivoire": "Cote_dIvoire",
    "Côte d'Ivoire": "Cote_dIvoire",
    "Ivory Coast": "Cote_dIvoire",
    "DRC": "Democratic_Republic_of_the_Congo",
    "Democratic Republic of Congo": "Democratic_Republic_of_the_Congo",
    "Democratic Republic of the Congo": "Democratic_Republic_of_the_Congo",
    "Congo Kinshasa": "Democratic_Republic_of_the_Congo",
    "Republic of the Congo": "Republic_of_the_Congo",
    "Congo Brazzaville": "Republic_of_the_Congo",
}

def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)

def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")

def clean_header_name(x: str) -> str:
    """Normalize TSV header names while preserving readable output names elsewhere."""
    return x.strip().replace("\ufeff", "")

def normalize_value(x: object, replace_space: bool = True) -> str:
    """Normalize metadata values for Nextstrain/Auspice categorical fields."""
    s = "" if x is None else str(x).strip()
    s = s.replace("\ufeff", "")
    if s.lower() in UNKNOWN_VALUES:
        return "unknown"
    s = s.replace("'", "").replace("’", "")
    s = s.replace("/", "_")
    s = s.replace("(", "").replace(")", "")
    s = s.replace("[", "").replace("]", "")
    s = s.replace("{", "").replace("}", "")
    s = s.replace("#", "")
    s = s.replace(">", "").replace("<", "")
    s = re.sub(r"\s+", "_" if replace_space else " ", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_") or "unknown"

def normalize_country(x: object) -> str:
    raw = "" if x is None else str(x).strip()
    if raw.lower() in UNKNOWN_VALUES:
        return "unknown"
    if raw in COUNTRY_ALIASES:
        return COUNTRY_ALIASES[raw]
    cleaned = normalize_value(raw)
    return COUNTRY_ALIASES.get(cleaned, cleaned)

def clean_date(x: object) -> str:
    """Convert GISAID dates to Augur-friendly uncertain dates."""
    d = "" if x is None else str(x).strip()
    d = d.replace("/", "-")
    if d.lower() in UNKNOWN_VALUES:
        return "XXXX-XX-XX"

    # Keep only the leading date-like part if GISAID adds notes.
    m = re.search(r"(\d{4})(?:-(\d{1,2}|XX))?(?:-(\d{1,2}|XX))?", d)
    if not m:
        return d

    year = m.group(1)
    month = m.group(2)
    day = m.group(3)

    if not month:
        return f"{year}-XX-XX"
    if month == "00":
        month = "XX"
    elif month != "XX":
        month = month.zfill(2)

    if not day:
        return f"{year}-{month}-XX"
    if day == "00":
        day = "XX"
    elif day != "XX":
        day = day.zfill(2)

    return f"{year}-{month}-{day}"

def split_gisaid_location(x: object) -> Tuple[str, str, str, str]:
    """Split strings like 'Africa / Madagascar / ...' into 4 geographic levels."""
    s = "" if x is None else str(x).strip()
    parts = [p.strip() for p in s.split("/") if p.strip()]
    while len(parts) < 4:
        parts.append("unknown")
    region = normalize_value(parts[0])
    country = normalize_country(parts[1])
    division = normalize_value(parts[2])
    location = normalize_value(parts[3])
    return region, country, division, location

def read_tsv_dicts(path: Path) -> Iterable[Dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return
        reader.fieldnames = [clean_header_name(x) for x in reader.fieldnames]
        for row in reader:
            yield {
                clean_header_name(k): ("" if v is None else str(v).strip())
                for k, v in row.items()
                if k is not None
            }

def get_field(row: Dict[str, str], *names: str, default: str = "") -> str:
    for name in names:
        if name in row:
            return row.get(name, default)
    # Fallback: case-insensitive lookup.
    low = {k.lower(): k for k in row}
    for name in names:
        key = low.get(name.lower())
        if key is not None:
            return row.get(key, default)
    return default

def parse_accession_from_header(header: str) -> Optional[str]:
    m = ACC_RE.search(header)
    if m:
        return m.group(0)
    # Fallback for 'virus|accession|date' headers if accession does not match EPI_ISL_ pattern.
    parts = header.split("|")
    if len(parts) >= 2 and parts[1].strip():
        return parts[1].strip()
    return None

def iter_fasta_records(path: Path) -> Iterable[Tuple[str, List[str]]]:
    header: Optional[str] = None
    seq: List[str] = []
    with open_text(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.upper())
        if header is not None:
            yield header, seq

def write_fasta_record(handle, name: str, seq_lines: Sequence[str], width: int = 80) -> None:
    seq = "".join(seq_lines).upper().replace(" ", "")
    handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")

def load_seqtech(seq_tech_path: Path) -> Dict[str, str]:
    seqtech: Dict[str, str] = {}
    for row in read_tsv_dicts(seq_tech_path):
        acc = get_field(row, "Accession ID", "accession")
        if not acc:
            continue
        val = get_field(row, "Sequencing technology", "seq_tech", default="unknown")
        seqtech[acc] = normalize_value(val) if val else "unknown"
    return seqtech

def build_patient_rows(patient_status_path: Path, seqtech: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    rows_by_acc: Dict[str, Dict[str, str]] = {}
    duplicates: List[str] = []

    for row in read_tsv_dicts(patient_status_path):
        acc = get_field(row, "Accession ID", "accession")
        virus_name = get_field(row, "Virus name", "virus_name")
        if not acc:
            continue
        if acc in rows_by_acc:
            duplicates.append(acc)
            continue

        region, country, division, location = split_gisaid_location(get_field(row, "Location", "location"))
        additional_location = normalize_value(get_field(row, "Additional location information", default="unknown"))

        rows_by_acc[acc] = {
            "strain": acc,
            "virus": "hChikV",
            "accession": acc,
            "virus_name": virus_name,
            "date": clean_date(get_field(row, "Collection date", "date")),
            "submission_date": clean_date(get_field(row, "Submission date", default="")),
            "region": region,
            "country": country,
            "division": division,
            "location": location,
            "additional_location_information": additional_location,
            "genotype": normalize_value(get_field(row, "Genotype", default="unknown")),
            "host": normalize_value(get_field(row, "Host", default="unknown")),
            "sex": normalize_value(get_field(row, "Gender", "sex", default="unknown")),
            "patient_age": normalize_value(get_field(row, "Patient age", default="unknown")),
            "patient_status": normalize_value(get_field(row, "Patient status", default="unknown")),
            "sampling_strategy": normalize_value(get_field(row, "Sampling strategy", default="unknown")),
            "passage": normalize_value(get_field(row, "Passage", default="unknown")),
            "specimen": normalize_value(get_field(row, "Specimen", default="unknown")),
            "seq_tech": seqtech.get(acc, "unknown"),
            "length": normalize_value(get_field(row, "Sequence Length", "length", default="unknown")),
        }

    if duplicates:
        eprint(f"Warning: duplicated accession IDs in patient_status.tsv; kept first occurrence: {len(duplicates)}")
    return rows_by_acc

def deterministic_hex(value: str) -> str:
    """Generate a stable readable-ish hex color from a string without using random."""
    digest = hashlib.md5(value.encode("utf-8")).hexdigest()
    # Avoid very dark colors by mixing each channel with 0x60.
    r = (int(digest[0:2], 16) + 0x60) // 2
    g = (int(digest[2:4], 16) + 0x60) // 2
    b = (int(digest[4:6], 16) + 0x60) // 2
    return f"#{r:02X}{g:02X}{b:02X}"

def write_colors(metadata_rows: List[Dict[str, str]], out_path: Path) -> None:
    countries = sorted({r["country"] for r in metadata_rows if r.get("country") and r.get("country") != "unknown"})
    genotypes = sorted({r["genotype"] for r in metadata_rows if r.get("genotype") and r.get("genotype") != "unknown"})
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for country in countries:
            writer.writerow(["country", country, deterministic_hex("country:" + country)])
        for genotype in genotypes:
            writer.writerow(["genotype", genotype, deterministic_hex("genotype:" + genotype)])

def is_unknown_geo_or_date(row: Dict[str, str]) -> bool:
    return (
        row.get("date", "").startswith("XXXX")
        or row.get("region", "unknown") == "unknown"
        or row.get("country", "unknown") == "unknown"
    )

def write_vip_color(metadata_rows: List[Dict[str, str]], out_path: Path) -> None:
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in metadata_rows:
            text = " ".join([
                row.get("country", ""),
                row.get("division", ""),
                row.get("location", ""),
                row.get("additional_location_information", ""),
                row.get("virus_name", ""),
            ]).lower()
            if "shenzhen" in text or re.search(r"\bsz\b", text):
                writer.writerow([row["strain"], "red"])
            elif row.get("country") == "China":
                writer.writerow([row["strain"], "blue"])

def command_prepare(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seqtech = load_seqtech(Path(args.seq_tech))
    rows_by_acc = build_patient_rows(Path(args.patient_status), seqtech)

    fasta_without_metadata: List[str] = []
    duplicated_fasta: List[str] = []
    written_accessions: List[str] = []
    seen: set[str] = set()

    sequences_path = outdir / "sequences.fasta"
    with open(sequences_path, "w", encoding="utf-8", newline="") as out_fasta:
        for header, seq_lines in iter_fasta_records(Path(args.fasta)):
            acc = parse_accession_from_header(header)
            if not acc:
                fasta_without_metadata.append(header)
                continue
            if acc not in rows_by_acc:
                fasta_without_metadata.append(acc)
                continue
            if acc in seen:
                duplicated_fasta.append(acc)
                continue
            seen.add(acc)
            written_accessions.append(acc)
            write_fasta_record(out_fasta, acc, seq_lines)

    metadata_rows = [rows_by_acc[acc] for acc in written_accessions]
    metadata_without_sequence = sorted(set(rows_by_acc) - set(written_accessions))

    metadata_cols = [
        "strain", "virus", "accession", "virus_name", "date", "submission_date",
        "region", "country", "division", "location", "additional_location_information",
        "genotype", "host", "sex", "patient_age", "patient_status",
        "sampling_strategy", "passage", "specimen", "seq_tech", "length",
    ]

    metadata_path = outdir / "metadata.tsv"
    with open(metadata_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in metadata_rows:
            writer.writerow(row)

    with open(outdir / "unknown.exclude", "w", encoding="utf-8", newline="") as handle:
        for row in metadata_rows:
            if is_unknown_geo_or_date(row):
                handle.write(row["strain"] + "\n")

    write_vip_color(metadata_rows, outdir / "vip.color.tsv")
    write_colors(metadata_rows, outdir / "colors.tsv")

    for name, values in {
        "fasta_without_metadata.txt": fasta_without_metadata,
        "duplicated_fasta_accessions.txt": duplicated_fasta,
        "metadata_without_sequence.txt": metadata_without_sequence,
    }.items():
        with open(outdir / name, "w", encoding="utf-8") as handle:
            for value in values:
                handle.write(str(value) + "\n")

    print(f"sequences.fasta: {len(written_accessions)} records")
    print(f"metadata.tsv: {len(metadata_rows)} rows")
    print(f"unknown.exclude: {sum(1 for r in metadata_rows if is_unknown_geo_or_date(r))} strains")
    print(f"fasta_without_metadata: {len(fasta_without_metadata)}")
    print(f"metadata_without_sequence: {len(metadata_without_sequence)}")
    print(f"duplicated_fasta_accessions: {len(duplicated_fasta)}")

def load_centroids(path: Path) -> Dict[Tuple[str, str], Tuple[str, str]]:
    centroids: Dict[Tuple[str, str], Tuple[str, str]] = {}
    with open_text(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 4:
                continue
            res, name, lat, lon = [x.strip() for x in row[:4]]
            if not res or not name or res.lower() in {"resolution", "geo_resolution"}:
                continue
            # Skip a possible header.
            if res.lower() == "country" and name.lower() in {"country", "name", "value"}:
                continue
            norm_name = normalize_country(name) if res == "country" else normalize_value(name)
            centroids[(res, norm_name)] = (lat, lon)
            # Also keep raw name for exact matching.
            centroids[(res, name)] = (lat, lon)
    return centroids

def command_add_coords(args: argparse.Namespace) -> None:
    metadata_path = Path(args.metadata)
    centroids_path = Path(args.centroids)
    out_path = Path(args.output)
    lat_longs_path = Path(args.lat_longs)
    missing_path = Path(args.missing_output)

    centroids = load_centroids(centroids_path)

    with open_text(metadata_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Empty metadata file: {metadata_path}")
        fieldnames = [clean_header_name(x) for x in reader.fieldnames]
        rows = [{clean_header_name(k): ("" if v is None else v.strip()) for k, v in row.items()} for row in reader]

    if "latitude" not in fieldnames:
        fieldnames.append("latitude")
    if "longitude" not in fieldnames:
        fieldnames.append("longitude")

    used_countries: Dict[str, Tuple[str, str]] = {}
    missing_countries: set[str] = set()

    for row in rows:
        country = normalize_country(row.get("country", "unknown"))
        row["country"] = country
        lat_lon = centroids.get(("country", country))
        if lat_lon is None and country != "unknown":
            missing_countries.add(country)
            row["latitude"] = ""
            row["longitude"] = ""
        elif lat_lon is not None:
            lat, lon = lat_lon
            row["latitude"] = lat
            row["longitude"] = lon
            used_countries[country] = (lat, lon)
        else:
            row["latitude"] = ""
            row["longitude"] = ""

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    # Augur lat-longs format: resolution, value, latitude, longitude. No header.
    with open(lat_longs_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for country in sorted(used_countries):
            lat, lon = used_countries[country]
            writer.writerow(["country", country, lat, lon])

    with open(missing_path, "w", encoding="utf-8") as handle:
        for country in sorted(missing_countries):
            handle.write(country + "\n")

    print(f"metadata with coords: {out_path}")
    print(f"lat_longs.tsv: {len(used_countries)} countries")
    print(f"missing country coordinates: {len(missing_countries)}")
    if missing_countries:
        eprint("Warning: some countries were not found in countries.centroids.tsv; see", missing_path)

def command_write_config(args: argparse.Namespace) -> None:
    config = {
        "title": args.title,
        "maintainers": [{"name": args.maintainer}],
        "colorings": [
            {"key": "country", "title": "Country", "type": "categorical"},
            {"key": "division", "title": "Division", "type": "categorical"},
            {"key": "region", "title": "Region", "type": "categorical"},
            {"key": "genotype", "title": "Genotype", "type": "categorical"},
            {"key": "host", "title": "Host", "type": "categorical"},
            {"key": "seq_tech", "title": "Sequencing technology", "type": "categorical"},
            {"key": "date", "title": "Date", "type": "continuous"},
        ],
        "filters": ["region", "country", "division", "genotype", "host", "seq_tech"],
        "display_defaults": {
            "color_by": "country",
            "geo_resolution": "country",
        },
        "panels": ["tree", "map", "entropy"],
    }

    out = Path(args.output)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"auspice config: {out}")

def command_count_values(args: argparse.Namespace) -> None:
    values = set()
    with open_text(Path(args.metadata), "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or args.column not in reader.fieldnames:
            print(0)
            return
        for row in reader:
            val = (row.get(args.column) or "").strip()
            if args.ignore_unknown and val.lower() in UNKNOWN_VALUES:
                continue
            values.add(val)
    print(len(values))

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare hChikV GISAID data for Nextstrain/Augur.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prepare", help="Create sequences.fasta, metadata.tsv, colors.tsv, and exclude files.")
    p.add_argument("--fasta", required=True, help="GISAID FASTA, e.g. gisaid/gisaid.fasta")
    p.add_argument("--patient-status", required=True, help="GISAID patient_status.tsv")
    p.add_argument("--seq-tech", required=True, help="GISAID seq_tech.tsv")
    p.add_argument("--outdir", required=True, help="Output working directory")
    p.set_defaults(func=command_prepare)

    p = sub.add_parser("add-coords", help="Add country centroids to filtered metadata and write lat_longs.tsv.")
    p.add_argument("--metadata", required=True, help="Input metadata, usually filtered_metadata.tsv")
    p.add_argument("--centroids", required=True, help="countries.centroids.tsv")
    p.add_argument("--output", default="filtered_metadata_with_coords.tsv")
    p.add_argument("--lat-longs", default="lat_longs.tsv")
    p.add_argument("--missing-output", default="missing_country_coords.tsv")
    p.set_defaults(func=command_add_coords)

    p = sub.add_parser("write-config", help="Write auspice_config.json once, outside the shell script.")
    p.add_argument("--output", default="auspice_config.json")
    p.add_argument("--title", default="hChikV genomic analysis")
    p.add_argument("--maintainer", default="Jie Huang")
    p.set_defaults(func=command_write_config)

    p = sub.add_parser("count-values", help="Count distinct values in a metadata column.")
    p.add_argument("--metadata", required=True)
    p.add_argument("--column", required=True)
    p.add_argument("--ignore-unknown", action="store_true")
    p.set_defaults(func=command_count_values)

    return parser

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
