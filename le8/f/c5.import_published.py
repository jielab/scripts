#!/usr/bin/env python3
"""Import Chen et al. Nature 2026 CIGMA Tables S3/S13, preserving each test family.

The original XLSX also retains Table S4 estimates. S4 has no inferential P values
and is not substituted for either cis-only test family. No models are fitted.
"""
import argparse
from collections import Counter
import gzip
import importlib.util
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("c5", Path(__file__).with_name("c5.cellulation.py"))
c5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c5)

ARTICLE = "https://www.nature.com/articles/s41586-026-10577-6"
TABLES = [
    ("Table S3 Onek1k cis only", "OneK1K", 10288, ["BIN", "BMem", "CD4ET", "CD4NC", "CD8ET", "CD8NC", "NK"]),
    ("Table S13 CLUES meta-analysis", "CLUES_ImmVar_meta", 10553, ["B", "NK", "T4", "T8", "cDC", "cM", "ncM"]),
]


def gene_map(annotation):
    mapping = {}
    opener = gzip.open if str(annotation).endswith(".gz") else open
    with opener(annotation, "rt") as stream:
        for line in stream:
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            if "gene_id" in attrs and "gene_name" in attrs:
                mapping[attrs["gene_id"].split(".")[0]] = c5.gene_name(attrs["gene_name"])
    counts = Counter(mapping.values())
    # Ambiguous symbols keep their stable Ensembl IDs; never collapse tests.
    return {key: value for key, value in mapping.items() if counts[value] == 1}


def run(workbook, annotation, output):
    output.mkdir(parents=True, exist_ok=True)
    names = gene_map(annotation)
    results, cells, records = [], [], []
    for sheet, cohort, size, cell_types in TABLES:
        original = pd.read_excel(workbook, sheet_name=sheet, header=2)
        if len(original) != size or original.gene.duplicated().any():
            raise ValueError(f"Unexpected complete family in {sheet}: {len(original)} rows")
        original.insert(0, "source_gene_id", original.gene)
        original["gene"] = original.source_gene_id.map(lambda gene: names.get(str(gene).split('.')[0], gene))
        if original.gene.duplicated().any():
            raise ValueError("Gene mapping would collapse separate source tests")
        # S13 reports cell estimates/SEs, but not individual-cell P values.
        original["cell_p_not_reported"] = np.nan
        directory = output / cohort
        directory.mkdir(exist_ok=True)
        source = directory / "published_table.csv"
        c5.write_csv(original, source)
        columns = dict(gene="gene", specific_p="p:V", shared_p="p:sigma_g2",
                       shared_variance="sigma_g2", shared_se="se:sigma_g2", specificity="specificity")
        if "v" in original:
            columns["mean_specific_variance"] = "v"
        cell_map = [dict(cell_type=cell, p_col=f"p:V_{cell}" if f"p:V_{cell}" in original else "cell_p_not_reported",
                         variance_col=f"V_{cell}", se_col=f"se:V_{cell}") for cell in cell_types]
        mapping = dict(source=ARTICLE, source_version="2026-08-26; " + sheet, cohort=cohort,
                       tissue="PBMC", component="cis", full_test_family=True, n_tests=size,
                       columns=columns, cells=cell_map, complete_cell_test_family=cohort == "OneK1K",
                       n_cell_tests=size * len(cell_types))
        mapfile = directory / "column_mapping.json"
        c5.write_json(mapping, mapfile)
        table = c5.import_cigma(source, mapfile, directory)
        table["source_gene_id"] = original.source_gene_id.to_numpy()
        # The paper defines cs-eGenes using Bonferroni, not the additional BH q.
        table["source_cs_egene"] = table.specific_p < .05 / size
        table["source_significance_rule"] = f"Published rule: P < 0.05/{size} (Bonferroni within {cohort})"
        c5.write_csv(table, directory / "cigma.results.csv")
        cell = c5.read_table(directory / "cigma.cell_types.csv")
        results.append(table)
        cells.append(cell)
        records.append(dict(sheet=sheet, cohort=cohort, genes=size, cell_types=cell_types,
                            source_cs_egenes=int(table.source_cs_egene.sum()),
                            mapped_symbols=int((original.gene != original.source_gene_id).sum()),
                            unmapped_ids=int((original.gene == original.source_gene_id).sum()),
                            gene_adjustment=f"BH across all {size} source genes; publication Bonferroni rule retained separately",
                            cell_adjustment="BH across complete gene x cell family" if cohort == "OneK1K" else "No cell P values supplied; none inferred"))
    c5.write_csv(pd.concat(results, ignore_index=True), output / "cigma.results.csv")
    c5.write_csv(pd.concat(cells, ignore_index=True), output / "cigma.cell_types.csv")
    provenance = dict(version=c5.VERSION, record_type="external_precomputed", source=ARTICLE,
                      source_version="2026-08-26", source_table=c5.stamp(workbook),
                      gene_annotation=c5.stamp(annotation), import_code=c5.stamp(Path(__file__)),
                      full_test_family=True, families=records,
                      gene_mapping="Unique Ensembl gene ID to GTF gene_name; unmapped/ambiguous IDs retained; no tests dropped",
                      table_s4="Retained in original workbook; estimates without inferential P values are not pooled with cis-only tests",
                      interpretation="Published PBMC expression-regulatory evidence; no UKB CIGMA fitting and no CAD-specific causal inference")
    c5.write_json(provenance, output / "cigma.provenance.json")
    print(json.dumps(records, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--gene-annotation", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    run(args.workbook, args.gene_annotation, args.outdir)
