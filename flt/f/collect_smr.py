"""Intersect same-tissue, same-probe cis-SMR results across the two traits."""
import csv
import sys
from pathlib import Path

run = Path(sys.argv[1])
traits = ["G732", "mdd_2025"]
fields = ["Tissue", "Gene_id", "Gene", "probeID"] + [f"{key}_{t}" for t in traits for key in ("p_SMR","p_HEIDI")]
with (run / "matched.tsv").open("w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    for tissue in ("eQTLGen","Whole_Blood"):
        kept = []
        for trait in traits:
            with (run / f"{tissue}_{trait}.smr").open() as src:
                keys = src.readline().split()
                if set(("probeID","Gene","p_SMR","p_HEIDI")) - set(keys):
                    raise ValueError(f"{tissue}/{trait}: invalid SMR header")
                passing = {}
                for line in src:
                    row = dict(zip(keys, line.split(), strict=True))
                    if row["p_SMR"] in ("NA","nan",".") or row["p_HEIDI"] in ("NA","nan","."):
                        continue
                    if float(row["p_SMR"]) < .05 and float(row["p_HEIDI"]) > .05:
                        probe = row["probeID"]
                        if probe in passing:
                            raise ValueError(f"Repeated probe {probe} in {tissue}/{trait}")
                        passing[probe] = row
                kept.append(passing)
        for probe in sorted(kept[0].keys() & kept[1].keys()):
            a,b = (result[probe] for result in kept)
            if a["Gene"] != b["Gene"]:
                raise ValueError(f"Inconsistent gene annotation at {probe}")
            record = dict(Tissue=tissue, Gene_id=probe.split(".")[0], Gene=a["Gene"], probeID=probe)
            for trait,row in zip(traits,(a,b)):
                record.update({f"{k}_{trait}": row[k] for k in ("p_SMR","p_HEIDI")})
            writer.writerow(record)
# Publish a complete result atomically. Run-specific originals remain available.
link = run / "publish.tmp"
link.symlink_to(run.name + "/matched.tsv")
link.replace(run.parent / "matched.tsv")
