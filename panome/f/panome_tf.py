#!/usr/bin/env python3
"""Fine-tune a pretrained numerical Transformer to predict disease Y from prot/met."""
import argparse
import copy
import json
import os
from pathlib import Path
import sys
from panome import parser as base_parser, configure
from tf_download import DEFAULT_MODEL, verify_model, download_model


def parser():
    base = base_parser()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Example: ./panome_TF.sh --Y cvd_cad,ra --biom prot,met")
    supported = """command trait biom ukb_phe phe_file omics_file input_source met_input met_map id_col
        baseline_col diagnosis_col death_col lost_col end_date disease_evidence_col healthy_date_cols
        covariates residualize categorical group_col split_file transform exclude_features feature_missing
        sample_missing horizon min_censor_survival min_events folds c_grid max_iter analysis_root run_name
        run_dir output seed cores r_bin max_samples demo_features demo train_only dry_run preflight full_input_hash
        replace device check_device batch_size epochs patience learning_rate weight_decay amp""".split()
    for action in base._actions:
        if action.dest in supported:
            action = copy.copy(action)
            if action.dest == "biom":
                action.choices = None
            p._add_action(action)
    p.set_defaults(analysis_root="/mnt/d/analysis/panome_TF", run_name="tabiclv2_finetuned",
        epochs=10, patience=3, learning_rate=1e-5, weight_decay=.01, batch_size=64)
    p.add_argument("--model-path", default=DEFAULT_MODEL)
    p.add_argument("--download-model", action="store_true", help="Download/verify the model and exit")
    p.add_argument("--max-features", type=int, default=256, help="Build-only feature screening; 0 retains all QC-passing assays")
    p.add_argument("--context-size", type=int, default=256, help="Maximum labeled build context rows")
    p.add_argument("--predict-batch-size", type=int, default=256)
    p.add_argument("--validation-samples", type=int, default=2048, help="Known tune rows used for early stopping; 0 uses all")
    p.add_argument("--with-clinical", action="store_true", help="Also supply baseline clinical covariates to the Transformer")
    p.add_argument("--unfreeze-encoder", action="store_true", help="Fine-tune column/row encoders as well as the ICL Transformer; uses more VRAM")
    return p


def parse_args(argv=None):
    # configure() shares the stable input and outcome rules of panome.sh;
    # its unrelated analysis defaults are internal, not exposed as TF options.
    p = parser()
    defaults = vars(base_parser().parse_args([]))
    defaults.update(p._defaults)
    return p.parse_args(argv, namespace=argparse.Namespace(**defaults))


def batch_configs(args):
    outcomes, layers = args.trait.split(","), args.biom.split(",")
    if any(not s.strip() or s != s.strip() for s in outcomes) or len(set(outcomes)) != len(outcomes):
        raise ValueError("--Y requires unique, nonempty, comma-separated outcomes")
    if not set(layers) <= {"prot", "met"} or len(set(layers)) != len(layers):
        raise ValueError("--biom expects prot, met, or prot,met")
    if len(layers) > 1 and args.omics_file:
        raise ValueError("--omics-file requires a single biom")
    if len(outcomes) > 1 and args.diagnosis_col:
        raise ValueError("--diagnosis-col requires a single Y")
    if len(outcomes)*len(layers) > 1 and (args.run_dir or args.output or args.command != "run"):
        raise ValueError("Batch mode requires separate default run directories and command run")
    for trait in outcomes:
        for biom in layers:
            a = copy.deepcopy(args)
            a.trait, a.biom = trait, biom
            a = configure(a)
            if a.cores < 1 or a.max_features < 0 or a.context_size < 4 or a.predict_batch_size < 1 or a.validation_samples < 0:
                raise ValueError("Invalid TF resource settings")
            if a.validation_samples and a.validation_samples < 2*a.min_events:
                raise ValueError("--validation-samples must be 0 or at least twice --min-events")
            yield a


def main(argv=None):
    args = parse_args(argv)
    configs = list(batch_configs(args))
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ[name] = str(args.cores)
    if args.download_model:
        print(json.dumps(download_model(args.model_path), indent=2))
        return
    if args.dry_run:
        for a in configs:
            out = Path(a.run_dir) if a.run_dir else Path(a.analysis_root)/a.trait/a.biom/a.run_name
            print(json.dumps(dict(output=str(out), model=a.model_path, config=vars(a)), indent=2))
        return
    import torch
    from neural import device_for
    from common import run_lock, prepare_run_directory, log, dump
    from threadpoolctl import threadpool_limits
    from tf_pipeline import train, evaluate, project
    if args.check_device:
        print(json.dumps(dict(torch=torch.__version__, cuda_available=torch.cuda.is_available(),
                             gpu=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None), indent=2))
        return
    import importlib.metadata
    if importlib.metadata.version("tabicl") != "2.2.0":
        raise RuntimeError("This adapter requires tabicl==2.2.0; install requirements-tf.txt using PANOME_PYTHON")
    for a in configs:
        a.device = device_for(a.device)
        out = Path(a.run_dir) if a.run_dir else Path(a.analysis_root)/a.trait/a.biom/a.run_name
        if a.command == "project":
            if not a.run_dir or not a.output or not args.phe_file or not args.omics_file:
                raise ValueError("project requires --run-dir, --phe-file, --omics-file and --output")
            with threadpool_limits(a.cores):
                project(a, out)
            continue
        if a.command == "evaluate":
            with run_lock(out), threadpool_limits(a.cores):
                evaluate(out)
            continue
        info = verify_model(a.model_path)
        paths = [] if a.demo else [a.phe_file, a.omics_file]
        paths += [v for v in [a.split_file, a.met_map if a.biom == "met" and a.met_input == "raw" and not a.demo else ""] if v]
        missing = [path for path in paths if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError("Missing inputs: "+", ".join(missing))
        if a.preflight:
            log("DONE", "TF preflight", f"Y={a.trait}; biom={a.biom}; device={a.device}; model checksum verified; data values not yet read")
            continue
        out.mkdir(parents=True, exist_ok=True)
        with run_lock(out), threadpool_limits(a.cores):
            if not prepare_run_directory(out, replace=a.replace, train_only=a.train_only):
                continue
            try:
                train(a, out, info)
                if not a.train_only:
                    evaluate(out)
            except Exception as exc:
                dump(out/"FAILED.json", dict(error_type=type(exc).__name__, message=str(exc)))
                raise


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[PANOME TF] ERROR {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        raise
