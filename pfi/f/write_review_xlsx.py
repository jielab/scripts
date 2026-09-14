#!/usr/bin/env python3
"""Write a simple Excel review workbook from the paper-score CSV.

This intentionally uses only the Python standard library so it works inside the
existing conda environment without adding openpyxl/xlsxwriter.
"""


# 🚩 Imports, paths, and shared inputs
from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

def col_name(n: int) -> str:
    out = ""
    while n:
        n, r = divmod(n - 1, 26)
        out = chr(65 + r) + out
    return out

def cell_ref(row: int, col: int) -> str:
    return f"{col_name(col)}{row}"

def is_number(value: str) -> bool:
    return bool(re.fullmatch(r"-?(?:\d+|\d+\.\d+)", value or ""))

def sheet_xml(rows: list[list[str]], sheet_name: str) -> str:
    max_cols = max((len(r) for r in rows), default=1)
    max_rows = max(len(rows), 1)
    dim = f"A1:{cell_ref(max_rows, max_cols)}"
    width = {
        1: 12, 2: 14, 3: 22, 4: 12, 5: 10, 6: 38, 7: 24, 8: 20,
        9: 12, 10: 12, 11: 12, 18: 12, 19: 30, 20: 80,
    }
    cols = "".join(
        f'<col min="{i}" max="{i}" width="{width.get(i, 12)}" customWidth="1"/>'
        for i in range(1, max_cols + 1)
    )
    body = []
    for r_idx, row in enumerate(rows, 1):
        cells = []
        for c_idx in range(1, max_cols + 1):
            value = row[c_idx - 1] if c_idx <= len(row) else ""
            ref = cell_ref(r_idx, c_idx)
            if r_idx > 1 and is_number(value):
                cells.append(f'<c r="{ref}"><v>{value}</v></c>')
            else:
                cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>')
        body.append(f'<row r="{r_idx}">{"".join(cells)}</row>')
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <dimension ref="{dim}"/>
  <sheetViews><sheetView workbookViewId="0"><pane xSplit="4" ySplit="1" topLeftCell="E2" activePane="bottomRight" state="frozen"/><selection pane="bottomRight" activeCell="E2" sqref="E2"/></sheetView></sheetViews>
  <cols>{cols}</cols>
  <sheetData>{"".join(body)}</sheetData>
</worksheet>'''

def write_xlsx(csv_path: Path, xlsx_path: Path, sheet_name: str) -> None:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = [row for row in csv.reader(f)]
    xlsx_path.parent.mkdir(parents=True, exist_ok=True)
    now = dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    files = {
        "[Content_Types].xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>''',
        "_rels/.rels": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>''',
        "docProps/core.xml": f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:creator>PFI pipeline</dc:creator>
  <cp:lastModifiedBy>PFI pipeline</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>
</cp:coreProperties>''',
        "docProps/app.xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>PFI pipeline</Application>
</Properties>''',
        "xl/workbook.xml": f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets><sheet name="{escape(sheet_name)}" sheetId="1" r:id="rId1"/></sheets>
</workbook>''',
        "xl/_rels/workbook.xml.rels": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>''',
        "xl/styles.xml": '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="1"><fill><patternFill patternType="none"/></fill></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
</styleSheet>''',
        "xl/worksheets/sheet1.xml": sheet_xml(rows, sheet_name),
    }
    with zipfile.ZipFile(xlsx_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet_name", default="TrainingData")
    args = ap.parse_args()
    write_xlsx(Path(args.csv), Path(args.xlsx), args.sheet_name)
    print(f"Wrote review xlsx: {args.xlsx}")

if __name__ == "__main__":
    main()
