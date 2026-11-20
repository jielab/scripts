#!/usr/bin/env python3
"""Export XLSX worksheet values as typed JSON for the R table packer.

Use the same ElementTree value-reading approach as 0.common.py, without
importing the analysis dispatcher.  In particular, inline text attributes
are XML metadata and boolean 0/1 values must remain distinct.  No model is
loaded, formula is executed, or third-party Python package is required.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import sys
import zipfile
from xml.etree import ElementTree as ET


NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL_ID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def workbook_values(path, max_bytes=None):
    """Read cached cell values; retain text, booleans, numbers and missing cells."""
    sheets = []
    with zipfile.ZipFile(path) as archive:
        if max_bytes is not None and sum(item.file_size for item in archive.infolist()) > max_bytes:
            raise ValueError("Workbook exceeds expanded size cap")
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(item.itertext()) for item in root.findall("s:si", NS)]
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {item.attrib["Id"]: item.attrib["Target"] for item in rels}
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        for sheet in book.findall("s:sheets/s:sheet", NS):
            target = targets[sheet.attrib[REL_ID]]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            if ".." in Path(target).parts:
                raise ValueError("Unsafe worksheet path")
            root = ET.fromstring(archive.read(target))
            rows = []
            for row in root.findall("s:sheetData/s:row", NS):
                values = []
                for cell in row.findall("s:c", NS):
                    reference = re.match(r"[A-Z]+", cell.attrib["r"])
                    if reference is None:
                        raise ValueError("Invalid worksheet cell reference")
                    column = 0
                    for letter in reference.group():
                        column = column * 26 + ord(letter) - 64
                    values.extend([None] * max(0, column - len(values)))
                    kind = cell.attrib.get("t")
                    value_node = cell.find("s:v", NS)
                    if kind == "inlineStr":
                        text_node = cell.find("s:is", NS)
                        value = "".join(text_node.itertext()) if text_node is not None else None
                    elif value_node is None or value_node.text is None:
                        value = None
                    elif kind == "s":
                        value = shared[int(value_node.text)]
                    elif kind == "b":
                        if value_node.text not in ("0", "1"):
                            raise ValueError("Invalid worksheet boolean")
                        value = value_node.text == "1"
                    elif kind in ("str", "e"):
                        value = value_node.text
                    else:
                        try:
                            value = float(value_node.text)
                            if not math.isfinite(value):
                                raise ValueError("Non-finite worksheet number")
                        except ValueError:
                            value = value_node.text
                    values[column - 1] = value
                if any(value is not None for value in values):
                    rows.append(values)
            columns = []
            if rows:
                width = max(map(len, rows))
                header = rows[0] + [None] * (width - len(rows[0]))
                body = [row + [None] * (width - len(row)) for row in rows[1:]]
                seen = {}
                for index, label in enumerate(header):
                    label = str(label) if label is not None else f"column_{index + 1}"
                    seen[label] = seen.get(label, 0) + 1
                    name = label if seen[label] == 1 else f"{label}__{seen[label]}"
                    values = [row[index] for row in body]
                    present = [value for value in values if value is not None]
                    if not present or all(isinstance(value, bool) for value in present):
                        kind = "logical"
                    elif all(isinstance(value, (int, float)) and not isinstance(value, bool) for value in present):
                        kind = "numeric"
                    else:
                        kind = "character"
                        values = [None if value is None else
                                  "TRUE" if value is True else "FALSE" if value is False else str(value)
                                  for value in values]
                    columns.append({"name": name, "type": kind, "values": values})
            sheets.append({"name": sheet.attrib["name"], "columns": columns})
    return {"format": "le8-xlsx-values-v1", "sheets": sheets}


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: xlsx_values.py INPUT.xlsx OUTPUT.json")
    result = workbook_values(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
