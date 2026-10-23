"""Reclaim completed GRID staging data; preview by default, --apply to release."""
import argparse
from contextlib import ExitStack
import fcntl
import gzip
import json
import subprocess
from pathlib import Path

from resume import atomic_json


def posterior_ready(run, root):
    """Bind historical chromosome moments to the published covariance table."""
    import importlib.util
    try:
        config = json.loads((run / "config.json").read_text())
        sources = [Path(x[0]) for x in config["files"] if Path(x[0]).name == "joint_posterior.h5"]
        traits = {path.relative_to(root / "csx").parts[0] for path in sources}
        if len(traits) != 1:
            return False
        provenance = json.loads((root / "csx" / traits.pop() / "provenance.json").read_text())
        output = Path(provenance["score_file"]).with_name("1csx.posterior.rds")
        if not output.is_file():
            return False
        spec = importlib.util.spec_from_file_location("grid_cleanup_io", Path(__file__).with_name("0.common.py"))
        io = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(io)
        metadata = io.rds_metadata(output)
        return metadata.get("schema") == "grid_csx_moments_v1" and json.loads(metadata["signature"]) == config
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError) as exc:
        print(f"[GRID cleanup] RETAIN posterior scratch {run}: {exc}", flush=True)
        return False


def merged_ready(path):
    """Require a complete, nonempty gzip stream before dropping its inputs."""
    try:
        with gzip.open(path, "rb") as stream:
            if not stream.readline() or not stream.readline():
                return False
            while stream.read(4 * 1024**2):
                pass
        return True
    except (OSError, EOFError):
        return False


def clean_csx(root, *, apply=False):
    """Keep merged scores/logs; release chromosome text and unpinned LD copies."""
    from abm_core.storage import busy_reason
    from dosage_storage import _safe_file
    root = Path(root)
    total = 0
    ld = root / "cache/csx-ld"
    for path in sorted(ld.glob("ld-*.hdf5")):
        try:
            _safe_file(path, ld)
            with path.with_suffix(".lock").open("a") as pin:
                fcntl.flock(pin, fcntl.LOCK_EX | fcntl.LOCK_NB)
                metadata = path.with_suffix(".json")
                saved = json.loads(metadata.read_text())
                source = Path(saved["source"])
                if not source.is_file() or source.resolve() == path.resolve():
                    continue
                stat = source.stat()
                if stat.st_size != saved["bytes"] or stat.st_mtime_ns != saved["mtime_ns"]:
                    continue
                if busy_reason([path]):
                    continue
                size = path.stat().st_size
                if apply:
                    path.unlink()
                    metadata.unlink()
                total += size
        except (OSError, ValueError, KeyError) as exc:
            print(f"[GRID cleanup] RETAIN {path}: {exc}", flush=True)
    # Trait locks also protect legacy score writers. Never unlink lock inodes.
    tables = root / "cache/tables"
    try:
        with ExitStack() as stack:
            for path in sorted(tables.glob("*/*/run.lock")):
                stream = stack.enter_context(path.open("a"))
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            for run in sorted(tables.glob("*/*/run*")):
                if not run.is_dir() or not (run / "config.json").is_file():
                    continue
                targets = []
                if run.parent.name == "scores":
                    for pop in ("AFR", "EAS", "EUR", "SAS"):
                        merged = run / (pop + ".tsv.gz")
                        if merged_ready(merged):
                            targets += list(run.glob(pop + ".chr*.sscore"))
                elif run.parent.name in ("meta", "auto"):
                    merged = run / "combined.gz"
                    if merged_ready(merged):
                        targets += list(run.glob("chr*.sscore"))
                elif run.parent.name == "posterior_scores" and any(run.glob("chr*/moments.h5")):
                    if posterior_ready(run, root):
                        targets += list(run.glob("chr*/moments.h5"))
                if not targets or busy_reason(targets):
                    continue
                for path in targets:
                    _safe_file(path, tables)
                size = sum(p.stat().st_size for p in targets)
                if apply:
                    atomic_json(run / "SCORES_RELEASED.json", dict(
                        files={str(p.relative_to(run)): p.stat().st_size for p in targets},
                        reason="chromosome scores consolidated; retained combined/published outputs and source inputs"))
                    for path in targets:
                        path.unlink()
                total += size
    except BlockingIOError:
        print("[GRID cleanup] RETAIN CSx score scratch: active trait lock", flush=True)
    print(f"[GRID cleanup] CSx {'released' if apply else 'candidates'}: {total / 2**30:.2f} GiB", flush=True)
    return total


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--config", default=str(Path(__file__).resolve().parents[1] / "config/abm.blocks.json"))
    p.add_argument("--traits", default="height,ldl,t2dm")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)
    if args.apply and args.dry_run:
        p.error("--apply and --dry-run are mutually exclusive")
    import block_inputs as b
    from abm_core.grid_inputs import locate_blocks
    from abm_core.storage import clean_completed_support
    from dosage_storage import release_dosages, MARKER
    config = b.read_config(args.config)
    root = Path(args.analysis_root).resolve()
    if any(not Path(config[key]).resolve().is_relative_to(root) for key in ("cache_root", "prsformer_cache")):
        p.error("Block config cache paths must be inside --analysis-root")
    layout = b.resolve_block_files(config["ld_block_dir"], config["ld_block_files"])
    traits = args.traits.split(",")
    for trait in traits:
        if trait not in ("height", "ldl", "t2dm"):
            p.error("Unknown trait: " + trait)
        base, _, _, _, marker = b.cache_paths(config, layout, trait)
        if args.apply and marker.is_file() and not (base / MARKER).is_file():
            print(f"[GRID cleanup] Verify {trait}; retain 32 real dosage rows before release", flush=True)
            arr, *_ = locate_blocks(dict(block_config=args.config, execution={"device": "cpu"}), trait, allow_build=False)
            arr._mmap.close()
    total = release_dosages(args.config, traits=traits, analysis_root=args.analysis_root, apply=args.apply)
    total += clean_completed_support(args.analysis_root, apply=args.apply)
    total += clean_csx(args.analysis_root, apply=args.apply)
    print(f"[GRID cleanup] {'Released' if args.apply else 'Releasable with existing probes'}: {total / 2**30:.2f} GiB", flush=True)


if __name__ == "__main__":
    main()
