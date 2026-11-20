"""Select UKB field columns from the raw RAP export when a prepared extract is absent."""
import csv
import gzip
import os
import re
import sys
import tempfile
from pathlib import Path

name, output = sys.argv[1:]
ukb = Path(os.environ["FLT_UKB"])
target = Path(output)
fields = {line.split()[0] for line in (ukb / "common" / (name + ".lst")).read_text(encoding="utf-8-sig").splitlines() if line.strip() and not line.startswith("#")}
fields = {re.sub(r"^p?(\d+)(?:_.*)?$", r"p\1", key) for key in fields}
if name == "flt":
    fields.add("p20750")
raw = Path(os.getenv("FLT_RAW", str(ukb / "rap/raw/pheno.tab.gz")))
with gzip.open(raw, "rt", newline="") as source, tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as tmp:
    tmp.close()
    try:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        def selected(column):
            return column == "eid" or (column.split("_")[0] in fields and
                                      ("_i" not in column or re.search(r"_i0(?:_|$)", column)))
        idx = [i for i,column in enumerate(header) if selected(column)]

        if "eid" not in header or len(idx) < 2:
            raise ValueError(f"No usable eid and {name} fields in {raw}")
        with gzip.open(tmp.name, "wt", newline="") as dest:
            w = csv.writer(dest, delimiter="\t", lineterminator="\n")
            w.writerow([header[i] for i in idx])
            for row in reader:
                w.writerow([row[i] for i in idx])
        os.replace(tmp.name, target)
    finally:
        Path(tmp.name).unlink(missing_ok=True)
