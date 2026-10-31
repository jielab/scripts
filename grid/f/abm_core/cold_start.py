"""Prepare missing inputs from raw sources without training any predictor."""
from pathlib import Path
import json
import sys
import argparse

ROOT = Path(__file__).resolve().parents[2]


def benchmark_parser(name, description):
    p = argparse.ArgumentParser(prog="./4.prsformer.sh " + name, description=description)
    p.add_argument("--traits", default="height,ldl,t2dm",
                   help="Comma-separated traits; train independent models sequentially in this order")
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--config", default=str(ROOT / "config/abm.json"))
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--threads", type=int, default=2 if name == "prsformer" else 4)
    p.add_argument("--dry-run", action="store_true")
    return p


def benchmark_config(args):
    from .pipeline import parser, configuration
    options = parser().parse_args(["--config", args.config, "--traits", args.traits,
                                   "--device", args.device, "--threads", str(args.threads)])
    config = configuration(options)
    args.config, args.traits = options.config, options.traits
    return config


def benchmark_device(device):
    from .support import _backend
    selected, _ = _backend(device)
    if selected != "cuda":
        raise ValueError("PRSformer neighborhood attention requires CUDA; "
                         "use abm --device cpu for ABM on CPU")
    from prsformer import neighborhood_runtime
    neighborhood_runtime()
    return selected


def preparation_command(args, config):
    """Generate or reuse the shared roster inside the analysis directory."""
    block_config = json.loads(Path(config["block_config"]).read_text())
    analysis = Path(args.analysis_root)
    return ["bash", str(ROOT / "f/prsformer.entry.sh"), "--python", sys.executable,
            "--traits", "height,ldl,t2dm", "--prsformer-mode", "prepare", "--no-replace",
            "--cache-dir", str(analysis / "cache/abm"), "--out-root", str(analysis / "abm"),
            "--prsformer-root", str(analysis), "--prsformer-cache-dir", block_config["prsformer_cache"],
            "--threads", str(args.threads), "--prsformer-threads", str(args.threads)]


def ensure_inputs(args, config):
    prepared = Path(getattr(args, "prepared", None) or
                    Path(args.analysis_root) / "cache/abm/prepared.joblib")
    needs_blocks = (getattr(args, "predictor", None) == "blocknet" or
                    getattr(args, "representation", None) == "csx_ldblocks_fixed")
    # An exported side-feature table is sufficient for score-only ABM and
    # benchmark resume; individual benchmark adapters restore dosage if needed.
    if prepared.is_file() and not needs_blocks:
        return
    block_config = json.loads(Path(config["block_config"]).read_text())
    cache = Path(block_config["prsformer_cache"])
    required = ("prepare.json", "data.tsv.gz", "variants.tsv.gz", "genotypes.npy")
    if all((cache / name).is_file() for name in required):
        return
    if needs_blocks:
        import block_inputs as b
        bc = b.read_config(config["block_config"])
        layout = b.resolve_block_files(getattr(args, "ld_block_dir", None) or bc["ld_block_dir"],
                                       getattr(args, "ld_block_files", None) or bc["ld_block_files"])
        if all(b.cache_paths(bc, layout, trait)[-1].is_file() for trait in args.traits.split(",")):
            # Reuse the reviewed verifier; no block building or training here.
            from .grid_inputs import locate_blocks
            for trait in args.traits.split(","):
                array, *_ = locate_blocks(config, trait, getattr(args, "ld_block_files", None),
                                          getattr(args, "ld_block_dir", None), allow_build=False)
                array._mmap.close()
            return
    if getattr(args, "smoke", False):
        raise FileNotFoundError("SMOKE cannot prepare full inputs; missing " + str(prepared))
    raise FileNotFoundError(
        "GRID shared preparation is missing. Run:\n"
        "  ./4.prsformer.sh prepare --traits height,ldl,t2dm\n"
        "This prepares inputs only; it does not train PRSformer.\n"
        f"Expected shared cache: {cache}; prepared features: {prepared}. "
        "For a custom analysis root, also pass matching --analysis-root and --config.")
