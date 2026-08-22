#!/usr/bin/env python3
"""Normalize sample metadata and create population/sex sample lists.

The input must have a header. Column names are auto-detected, or can be
specified explicitly. Output is a tab-delimited table with columns:
sample, pop, super_pop, sex.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ALIASES = {
    "sample": ["sample", "sample_id", "sampleid", "iid", "id", "eid"],
    "pop": ["pop", "population", "population_code", "group", "ancestry_group"],
    "super_pop": ["super_pop", "superpopulation", "super_population", "ancestry", "race", "ethnicity"],
    "sex": ["sex", "gender", "genetic_sex", "biological_sex"],
}


def norm_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", x.strip().lower()).strip("_")


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    lines = [x for x in path.read_text(encoding="utf-8-sig").splitlines() if x.strip()]
    if not lines:
        raise SystemExit(f"Empty sample panel: {path}")
    delimiter = "\t" if "\t" in lines[0] else None
    if delimiter:
        reader = csv.DictReader(lines, delimiter="\t")
        rows = [{str(k): (v or "").strip() for k, v in r.items()} for r in reader]
        return list(reader.fieldnames or []), rows
    header = re.split(r"\s+", lines[0].strip())
    rows = []
    for line_no, line in enumerate(lines[1:], start=2):
        vals = re.split(r"\s+", line.strip())
        if len(vals) != len(header):
            raise SystemExit(f"Malformed whitespace table {path}:{line_no}: expected {len(header)} fields, found {len(vals)}")
        rows.append(dict(zip(header, vals)))
    return header, rows


def choose(headers: list[str], explicit: str | None, role: str, required: bool) -> str | None:
    lookup = {norm_name(h): h for h in headers}
    if explicit:
        key = norm_name(explicit)
        if key not in lookup:
            raise SystemExit(f"Requested {role} column {explicit!r} is absent. Columns: {headers}")
        return lookup[key]
    for alias in ALIASES[role]:
        if alias in lookup:
            return lookup[alias]
    if required:
        raise SystemExit(f"Cannot detect required {role} column. Columns: {headers}")
    return None


def normalize_sex(value: str) -> str:
    x = value.strip().lower()
    if x in {"male", "m", "1", "xy"}:
        return "male"
    if x in {"female", "f", "2", "xx"}:
        return "female"
    if x in {"", "na", "nan", "unknown", "0", ".", "-9"}:
        return "unknown"
    raise ValueError(value)


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    if not value or value in {".", ".."}:
        raise ValueError("empty/unsafe group name")
    return value


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--list-dir", type=Path)
    p.add_argument("--sample-col")
    p.add_argument("--pop-col")
    p.add_argument("--super-pop-col")
    p.add_argument("--sex-col")
    p.add_argument("--require-sex", action="store_true")
    p.add_argument("--allow-missing-pop", action="store_true")
    args = p.parse_args()

    headers, rows = read_table(args.input)
    sample_col = choose(headers, args.sample_col, "sample", True)
    pop_col = choose(headers, args.pop_col, "pop", not args.allow_missing_pop)
    super_col = choose(headers, args.super_pop_col, "super_pop", False)
    sex_col = choose(headers, args.sex_col, "sex", args.require_sex)

    out_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for i, row in enumerate(rows, start=2):
        sample = row[sample_col].strip()
        if not sample or sample.startswith("#"):
            continue
        if sample in seen:
            raise SystemExit(f"Duplicate sample ID {sample!r} in {args.input}")
        seen.add(sample)
        pop = row[pop_col].strip() if pop_col else "ALL"
        super_pop = row[super_col].strip() if super_col else pop
        if not pop:
            if args.allow_missing_pop:
                pop = super_pop or "ALL"
            else:
                raise SystemExit(f"Missing population for sample {sample}")
        if not super_pop:
            super_pop = pop
        pop = safe_name(pop)
        super_pop = safe_name(super_pop)
        try:
            sex = normalize_sex(row[sex_col]) if sex_col else "unknown"
        except ValueError:
            raise SystemExit(f"Unrecognized sex value {row[sex_col]!r} for sample {sample}")
        if args.require_sex and sex == "unknown":
            raise SystemExit(f"Sex is required but unknown for sample {sample}")
        out_rows.append({"sample": sample, "pop": pop, "super_pop": super_pop, "sex": sex})

    if not out_rows:
        raise SystemExit(f"No samples found in {args.input}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "pop", "super_pop", "sex"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    if args.list_dir:
        args.list_dir.mkdir(parents=True, exist_ok=True)
        groups: dict[str, list[str]] = {"ALL": []}
        sex_groups: dict[str, list[str]] = {"ALL.female": [], "ALL.male": []}
        pop_super: dict[str, str] = {"ALL": "ALL"}
        for row in out_rows:
            sample, pop, sex = row["sample"], safe_name(row["pop"]), row["sex"]
            groups["ALL"].append(sample)
            groups.setdefault(pop, []).append(sample)
            pop_super[pop] = row["super_pop"]
            if sex in {"female", "male"}:
                sex_groups[f"ALL.{sex}"].append(sample)
                sex_groups.setdefault(f"{pop}.{sex}", []).append(sample)
        for name, ids in {**groups, **sex_groups}.items():
            if not ids:
                continue
            (args.list_dir / f"{name}.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")
        with (args.list_dir / "sample_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
            handle.write("group\tn_samples\tsuper_pop\tsex\tsample_file\n")
            for name, ids in sorted({**groups, **sex_groups}.items()):
                if not ids:
                    continue
                base = name.rsplit(".", 1)[0] if name.endswith((".female", ".male")) else name
                sex = name.rsplit(".", 1)[1] if name.endswith((".female", ".male")) else "all"
                handle.write(f"{name}\t{len(ids)}\t{pop_super.get(base, 'ALL')}\t{sex}\t{args.list_dir / (name + '.txt')}\n")

    print(f"normalized_samples={len(out_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
