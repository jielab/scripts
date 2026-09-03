#!/usr/bin/env python3
"""Run g:Profiler enrichment with a custom assayed-protein background."""

import csv
import json
import sys
import urllib.request


def read_ids(path):
    with open(path, encoding="utf-8") as handle:
        return list(dict.fromkeys(x.strip() for x in handle if x.strip()))


def main():
    if len(sys.argv) != 4:
        raise SystemExit("usage: c1_enrich.py SIGNIFICANT_IDS BACKGROUND_IDS OUTPUT_CSV")
    query, background = read_ids(sys.argv[1]), read_ids(sys.argv[2])
    payload = {
        "organism": "hsapiens",
        "query": query,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        # Return ranked terms even when none passes FDR 0.05. The R plot marks
        # how many displayed terms are significant instead of producing a
        # visually blank figure for a biologically null result.
        "user_threshold": 1.0,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": background,
        "no_evidences": True,
    }
    request = urllib.request.Request(
        "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "LE8-C1/1.0"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        result = json.load(response).get("result", [])
    fields = [
        "source", "term_id", "term_name", "adjusted_p", "intersection_size",
        "term_size", "query_size",
    ]
    with open(sys.argv[3], "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in result:
            if row.get("p_value") is not None:
                writer.writerow({
                    "source": row.get("source", ""),
                    "term_id": row.get("native", ""),
                    "term_name": row.get("name", ""),
                    "adjusted_p": row.get("p_value", ""),
                    "intersection_size": row.get("intersection_size", 0),
                    "term_size": row.get("term_size", ""),
                    "query_size": row.get("query_size", len(query)),
                })


if __name__ == "__main__":
    main()
