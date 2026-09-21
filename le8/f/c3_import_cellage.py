#!/usr/bin/env python3
"""Import the authors' published CellAge marker mapping, preserving provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import pandas as pd


def build(repo, output):
    source = repo / "preprocessing/cell_type_mapping_update.csv"
    version = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    table = pd.read_csv(source)
    rows = []
    for _, row in table.iterrows():
        for gene in str(row["Marker Genes"]).split(","):
            gene = gene.strip().upper()
            if gene and gene != "NAN":
                rows.append(dict(gene=gene, cell_type=row["Original Cell types"],
                    source="Ding et al. CellAge author-published marker mapping; HPA-derived", source_version=version))
    out = pd.DataFrame(rows).drop_duplicates(["gene", "cell_type"])
    if out.empty or out.isna().any().any():
        raise ValueError("Empty or malformed CellAge mapping")
    output.mkdir(parents=True, exist_ok=True)
    out.to_csv(output/"atlas.csv", index=False)
    shutil.copy2(repo/"LICENSE", output/"LICENSE")
    provenance = dict(repository="https://github.com/dingdaisy/cellage", commit=version,
        source_file=str(source.relative_to(repo)), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        paper="https://doi.org/10.1038/s41591-026-04446-y",
        rule="Authors' supplied Marker Genes, not their SomaScan-only Somamer list",
        interpretation="Cell-enriched expression annotation, not secretion source or CIGMA inference",
        rows=len(out), genes=int(out.gene.nunique()), cell_types=int(out.cell_type.nunique()))
    (output/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    build(a.repo, a.outdir)
