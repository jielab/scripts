#!/usr/bin/env python3
"""Prepare the shared GRID/PRSformer cohort and genotype cache without training.

This command makes an existing implementation dependency explicit:
GRID may need the shared genotype matrix while constructing CSx-weighted LD
blocks, but GRID never needs PRSformer/PRSformer2 fitted models or predictions.
"""
from __future__ import annotations
import argparse
import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

from abm_core.pipeline import parser as abm_parser, configuration
from abm_core.cold_start import preparation_command


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--traits", default="height,ldl,t2dm")
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--config", default=str(Path(__file__).resolve().parents[1] / "config/abm.json"))
    p.add_argument("--device", choices=["auto","cpu","cuda"], default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    if a.threads < 1:
        p.error("--threads must be positive")

    # Reuse the reviewed ABM configuration parser so paths/traits stay identical.
    opts = abm_parser().parse_args([
        "--traits", a.traits,
        "--analysis-root", a.analysis_root,
        "--config", a.config,
        "--device", a.device,
        "--threads", str(a.threads),
    ])
    cfg = configuration(opts)
    cmd = preparation_command(opts, cfg)

    if a.dry_run:
        print("[4.prsformer.sh prepare] " + " ".join(map(str, cmd)))
        print("[4.prsformer.sh prepare] then verify/build cache/abm/prepared.joblib; no model training")
        return

    from kinship import analysis_guard
    with analysis_guard(a.analysis_root):
        env = dict(os.environ)
        env["GRID_CACHE_ROOT"] = str(Path(a.analysis_root) / "cache")
        env["PRSFORMER_KEEP_GENOTYPES"] = "1"
        env["GRID_KEEP_SHARED_INPUTS"] = "1"
        subprocess.run(cmd, check=True, env=env)

        # Create/verify GRID's small feature cache on exactly the same frozen roster.
        # This step does not train PRSformer or GRID predictors.
        from abm_prepare import prepare_abm
        prepare_abm(SimpleNamespace(
            analysis_root=a.analysis_root,
            cache_dir=str(Path(a.analysis_root) / "cache/abm"),
            prepared=None,
            traits="height,ldl,t2dm",
            out_root=str(Path(a.analysis_root) / "abm"),
            prepare_config=None,
        ))
        print("[4.prsformer.sh prepare] READY: shared inputs retained; no PRSformer model trained")



if __name__ == "__main__":
    main()
