#!/usr/bin/env python3
"""Manual real-VCF integration test; Atlas effects are explicitly simulated.

Run with the installed allele Python. Outputs never go into formal gnomAD results.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pandas as pd
from cyvcf2 import VCF

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="new test-only output directory")
    args = parser.parse_args()
    base = Path(args.out).resolve()
    base.mkdir(parents=True, exist_ok=False)
    inputs = base / "input_1kg_regions"
    inputs.mkdir()
    provenance = []
    for chrom in [1, 2, 22]:
        source = Path(f"/mnt/f/gen/1kg/38/vcf/chr{chrom}.vcf.gz")
        reader = VCF(str(source), lazy=True)
        samples = len(reader.samples)
        reader.close()
        region = f"chr{chrom}:20000001-21000000"
        target = inputs / f"chr{chrom}.vcf.gz"
        subprocess.run(
            [
                str(Path(sys.prefix) / "bin/bcftools"),
                "view",
                "--no-version",
                "--no-update",
                "--drop-genotypes",
                "--regions",
                region,
                "--output-type",
                "z",
                "--output-file",
                str(target),
                str(source),
            ],
            check=True,
        )
        provenance.append(
            dict(
                path=str(source),
                region=region,
                sample_count=samples,
                subset=str(target),
                counts="original cohort INFO retained; sample columns removed",
            )
        )
    output = base / "output_FAKE_ATLAS"
    (base / "NOTICE.txt").write_text(
        "Real 1KG genotypes summarized by original INFO. ALL AlphaGenome effects here are simulated for software testing, not biological results. No API used.\n"
    )
    with tempfile.TemporaryDirectory(prefix="allele-real-test-imports-") as tmp:
        inject = Path(tmp)
        (inject / "sitecustomize.py").write_text(
            f"import sys\nif sys.prefix == {sys.prefix!r}:\n    from alphagenome.atlas import atlas\n    from fake_atlas import FakeAtlas\n    atlas.create = lambda key: FakeAtlas()\n"
        )
        env = {
            **os.environ,
            "ALLELE_CONFIG": "/dev/null",
            "ALLELE_ENV_NAME": "allele",
            "ALPHA_GENOME_API_KEY": "offline-test-only",
            "PYTHONPATH": str(inject) + os.pathsep + str(ROOT / "tests"),
        }
        with (base / "pipeline.log").open("w") as log:
            subprocess.run(
                [
                    str(ROOT / "allele.sh"),
                    "all",
                    "--vcf-dir",
                    str(inputs),
                    "--fasta",
                    "/mnt/f/gen/1hgp/38/GRCH38.fasta.gz",
                    "--out",
                    str(output),
                    "--chunk",
                    "200",
                ],
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
    sites = pd.read_parquet(output / "work/sites.parquet")
    candidates = pd.read_parquet(output / "work/candidates.parquet")
    scores = pd.read_parquet(output / "results/stage1_candidates.parquet")
    assert len(candidates) == 3 * len(sites)
    assert (
        candidates.groupby("site_id")["observed"]
        .sum()
        .sort_index()
        .equals(sites.set_index("site_id").occupancy.sort_index())
    )
    assert (
        scores.loc[scores.occupancy < 3, "mutation_prediction_source"]
        .eq("out_of_fold")
        .all()
    )
    metrics = json.loads((output / "results/stage1_metrics.json").read_text())
    assert metrics["cv_group"] == "chrom" and metrics["cv_status"] == "ok"
    assert len(list((output / "final").glob("Fig*.png"))) == 4
    assert len(list((output / "final").glob("Table*.tsv"))) == 6
    summary = dict(
        inputs=provenance,
        site_count=len(sites),
        candidate_count=len(candidates),
        occupancy_counts={
            str(k): int(v)
            for k, v in sites.occupancy.value_counts().sort_index().items()
        },
        selected_sites=len(pd.read_parquet(output / "work/atlas_sites.parquet")),
        completed_stages="prepare select stage1 atlas geometry stage2 final".split(),
        stage1_metrics=metrics,
        atlas="SIMULATED; no real API; downstream effects/results are not biological evidence",
        output_directory=str(output),
    )
    (base / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
