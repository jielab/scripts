#!/usr/bin/env python3
"""Check imports and executables without contacting Atlas."""
import argparse
import importlib
import os
from pathlib import Path
import shutil
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    print("Python:", sys.executable)
    print("root:", args.root)
    print("out:", args.out)
    failed = []
    for command in ["bcftools", "samtools", "tabix"]:
        path = shutil.which(command)
        print(f'{command}: {path or "MISSING"}')
        if not path:
            failed.append(command)
    for name in [
        "numpy",
        "pandas",
        "scipy",
        "sklearn",
        "statsmodels",
        "statsmodels.discrete.conditional_models",
        "matplotlib",
        "pyarrow",
        "zarr",
        "numcodecs",
        "cyvcf2",
        "pysam",
        "tqdm",
        "alphagenome.atlas.atlas",
    ]:
        try:
            module = importlib.import_module(name)
            if name == "zarr" and int(module.__version__.split(".")[0]) != 2:
                raise RuntimeError("this pipeline requires zarr >=2,<3")
            print(f'{name}: OK {getattr(module, "__version__", "")}')
        except Exception as error:
            failed.append(name)
            print(f"{name}: FAILED ({error})")
    print(
        "ALPHA_GENOME_API_KEY:",
        (
            "SET"
            if os.getenv("ALPHA_GENOME_API_KEY")
            else "NOT SET (needed only for atlas)"
        ),
    )
    for name in ["ALLELE_VCF_DIR", "ALLELE_FASTA"]:
        value = os.getenv(name)
        present = bool(value and Path(value).exists())
        print(
            f'{name}: {value or "NOT SET"} [{"FOUND" if present else "MISSING; configure before prepare"}]'
        )
    if failed:
        print("Dependency check FAILED: " + ", ".join(failed), file=sys.stderr)
        return 1
    print("Dependencies OK. Data paths and API access must be configured separately.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
