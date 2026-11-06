"""Serial benchmark workflow. Publish only after every requested module succeeds."""
import json
import os
from pathlib import Path

from . import pipeline


def main(argv=None):
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = pipeline.parser()
    parser.prog = "./grid.sh workflow"
    parser.description = "Run prsformer, prsformer2, ABM, then final."
    args = parser.parse_args(argv)
    config = pipeline.configuration(args)
    abm_args = [arg for arg in argv if arg != "--module-only"]
    if args.module_only:
        return pipeline.main(abm_args)
    if args.smoke or args.through != "report" or args.output_root or args.report_from:
        parser.error("--smoke, --report-from, partial --through and --output-root belong to ./grid.sh abm")
    if args.dry_run:
        print(json.dumps(dict(
            modules=["prsformer", "prsformer2", "abm", "final"],
            before_modules=["self-check", "prepare_inputs_if_missing"],
            traits=args.traits.split(","), predictor=args.predictor,
            representation=args.representation, device=args.device,
            config_path=args.config, abm_phases=list(pipeline.PHASES),
            after_all_modules="final summary and result index",
        ), indent=2))
        return
    from resume import run_lock
    from .cold_start import ensure_inputs, benchmark_device
    from prsformer_resume import dispatch
    from .prsformer_adapter import main as prsformer2
    from final_dispatch import main as final
    benchmark_device(args.device)
    with run_lock(Path(args.analysis_root) / "workflow"):
        os.environ["PRSFORMER_KEEP_GENOTYPES"] = "1"
        os.environ.setdefault("GRID_CACHE_ROOT", str(Path(args.analysis_root) / "cache"))
        ensure_inputs(args, config)
        shared = ["--traits", args.traits, "--analysis-root", args.analysis_root,
                  "--device", args.device, "--threads", str(args.threads),
                  "--config", args.config]
        print("[GRID] MODULE 1/4 prsformer", flush=True)
        dispatch(pipeline.ROOT, shared)
        print("[GRID] MODULE 2/4 prsformer2", flush=True)
        prsformer2(shared)
        print("[GRID] MODULE 3/4 abm", flush=True)
        pipeline.main(abm_args)
        print("[GRID] MODULE 4/4 final: all requested modules completed", flush=True)
        final(["--analysis-root", args.analysis_root, "--traits", args.traits])
