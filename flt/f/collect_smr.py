"""Collect same-tissue, same-probe SMR evidence and retain all tested probes."""
import csv
import math
import sys
from pathlib import Path

TRAITS = ("G732", "mdd_2025")
FIELDS = ("Tissue", "Gene_id", "Gene", "probeID")


def pvalue(text):
    try:
        value = float(text)
    except (ValueError, TypeError):
        return None
    return value if math.isfinite(value) and 0 <= value <= 1 else None


def read_results(path):
    with path.open() as source:
        header = source.readline().split()
        missing = {"probeID", "Gene", "p_SMR", "p_HEIDI"} - set(header)
        if missing:
            raise ValueError(f"{path}: missing {', '.join(sorted(missing))}")
        rows = [dict(zip(header, line.split(), strict=True)) for line in source if line.strip()]
    if len({row["probeID"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: repeated probe ID")
    n = sum(pvalue(r["p_SMR"]) is not None for r in rows)
    for row in rows:
        p = pvalue(row["p_SMR"])
        row["p_Bonferroni"] = min(p * n, 1) if p is not None else None
        row["n_tested_SMR"] = n
        row["p_SMR"] = p
        row["p_HEIDI"] = pvalue(row["p_HEIDI"])
    return rows


def publish(path, fields, rows):
    temp = path.with_suffix(".part")
    with temp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)


def main(run, dest):
    tissues = (run / "tissues.txt").read_text().split()
    if not tissues or set(tissues) - {"eQTLGen", "Whole_Blood"}:
        raise ValueError("No valid completed SMR tissues")
    matched = []
    tests = []
    for tissue in tissues:
        candidates = []
        for trait in TRAITS:
            rows = read_results(run / f"{tissue}-{trait}.smr")
            tests.extend(dict(r, Tissue=tissue, Trait=trait) for r in rows)
            candidates.append({row["probeID"]: row for row in rows
                               if row["p_SMR"] is not None and row["p_SMR"] < .05
                               and row["p_HEIDI"] is not None and row["p_HEIDI"] > .05})
        for probe in sorted(candidates[0].keys() & candidates[1].keys()):
            a, b = [x[probe] for x in candidates]
            if a["Gene"] != b["Gene"]:
                raise ValueError(f"Inconsistent gene annotation at {probe}")
            row = dict(Tissue=tissue, Gene_id=probe.split(".")[0], Gene=a["Gene"], probeID=probe)
            for trait, fields in zip(TRAITS, [a, b]):
                row.update({f"{name}_{trait}": fields[name]
                            for name in ("p_SMR", "p_HEIDI", "p_Bonferroni")})
            row["both_bonferroni_pass"] = all(x["p_Bonferroni"] < .05 for x in (a, b))
            matched.append(row)
    all_fields = list(dict.fromkeys(["Tissue", "Trait", "probeID", "Gene", "p_SMR", "p_HEIDI",
                                     "p_Bonferroni", "n_tested_SMR"] + [k for row in tests for k in row]))
    publish(dest / "probe-tests.tsv", all_fields, tests)
    fields = list(FIELDS) + [f"{key}_{trait}" for trait in TRAITS
                            for key in ("p_SMR", "p_HEIDI", "p_Bonferroni")] + ["both_bonferroni_pass"]
    publish(dest / "matched.tsv", fields, matched)
    (dest / "tissues.txt").write_text("\n".join(tissues) + "\n")


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))
