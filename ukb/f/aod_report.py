#!/usr/bin/env python3
"""Render the AOD aggregate summary as the offline date-check.html report.

No individual-level input is read by this renderer. A small defensive gate rejects
obvious identifier/event payloads; aggregation and disclosure review belong to
the upstream analysis. This gate is not a statistical disclosure-control policy.
"""
from __future__ import annotations

import argparse
import html
import json
import os
from pathlib import Path
import tempfile


FORBIDDEN_KEYS = frozenset({
    "eid", "eids", "iid", "participant_id", "participant_ids", "patient_id",
    "patient_ids", "individual_records", "raw_events", "date_raw",
})


def validate_summary(data):
    if not isinstance(data, dict):
        raise ValueError("The report input must be one aggregate JSON object")
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported AOD report schema; expected schema_version=1")
    if not isinstance(data.get("trait"), str) or not data["trait"].strip():
        raise ValueError("The report input needs a non-empty trait")
    for key in ("cohort", "consistency", "burden", "clusters"):
        if not isinstance(data.get(key), dict):
            raise ValueError(f"The report input needs an aggregate {key} object")

    def visit(value, path="summary"):
        if isinstance(value, dict):
            for key, item in value.items():
                if key.lower() in FORBIDDEN_KEYS:
                    raise ValueError(f"Individual-level key is not permitted: {path}.{key}")
                visit(item, f"{path}.{key}")
        elif isinstance(value, list):
            for item in value:
                visit(item, path + "[]")

    visit(data)
    # Reject NaN/Infinity rather than silently converting an invalid statistic.
    json.dumps(data, ensure_ascii=False, allow_nan=False)


def render(data, output):
    validate_summary(data)
    template = Path(__file__).with_name("aod.html").read_text(encoding="utf-8")
    if template.count("__AOD_JSON__") != 1:
        raise ValueError("The aod.html template must contain exactly one aggregate JSON placeholder")
    packed = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    # JSON lives in a non-executable script element. Escape all HTML-sensitive
    # characters so a title or annotation cannot terminate that element.
    packed = (packed.replace("&", "\\u0026").replace("<", "\\u003c")
              .replace(">", "\\u003e").replace("\u2028", "\\u2028")
              .replace("\u2029", "\\u2029"))
    document = (template.replace("__TITLE__", html.escape(data["trait"], quote=True))
                .replace("__DEMO_HIDDEN__", "" if data.get("demo") else "hidden")
                .replace("__AOD_JSON__", packed))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent,
                                     prefix=".aod-report-", suffix=".tmp", delete=False) as handle:
        temp = Path(handle.name)
        try:
            handle.write(document)
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temp.unlink(missing_ok=True)
            raise
    try:
        os.replace(temp, output)
    finally:
        temp.unlink(missing_ok=True)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Aggregate summary.json from aod.py")
    parser.add_argument("--output", required=True, help="Self-contained date-check.html destination")
    args = parser.parse_args(argv)
    with open(args.input, encoding="utf-8") as stream:
        data = json.load(stream)
    result = render(data, args.output)
    print(f"AOD report: {result}")


if __name__ == "__main__":
    main()
