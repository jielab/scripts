"""Validate MRI sources before costly analyses or extraction begins."""
import csv
import gzip
from pathlib import Path


def header(path):
    path = Path(path)
    if not path.is_file() or not path.stat().st_size:
        raise FileNotFoundError(f"MRI input missing or empty: {path}. Prepare the YGB input first or set DOSE_MRI_RAW to an existing imaging table.")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        return handle.readline().rstrip("\r\n").split("\t")


def check_inputs(ygb, dictionary, raw=None, mode="all"):
    ygb = Path(ygb)
    raw = Path(raw) if raw else ygb / "rap/img.tab.gz"
    fields = []
    if mode in ("mri", "all"):
        available = set(header(raw))
        with open(dictionary, encoding="utf-8-sig", newline="") as handle:
            requested = list(dict.fromkeys("p" + row["field_id"] + "_i2" for row in csv.DictReader(handle)))
        fields = [field for field in requested if field in available]
        if not fields:
            raise ValueError(f"{raw}: no MRI dictionary fields present")
        missing = {"eid", "p26521_i2"} - available
        if missing:
            raise ValueError(f"{raw}: missing required MRI fields {sorted(missing)}")
        fields = list(dict.fromkeys(["eid", *fields, "p26521_i2"]))
    if mode in ("visit", "all"):
        visit = ygb / "rap/raw/pheno.tab.gz"
        missing = {"eid", "p21003_i2", "p54_i2", "p53_i2"} - set(header(visit))
        if missing:
            raise ValueError(f"{visit}: missing imaging visit fields {sorted(missing)}")
    return raw, fields
