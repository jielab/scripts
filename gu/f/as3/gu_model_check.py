#!/usr/bin/env python3
"""Load the bundled AS3 model definitions and external checkpoint pair."""
from __future__ import annotations

import argparse
import pickle
import sys
import warnings
from pathlib import Path


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.runtime.resolve()))
    warnings.filterwarnings(
        "ignore",
        message=r"enable_nested_tensor is True, but self\.use_nested_tensor is False.*",
        category=UserWarning,
        module=r"torch\.nn\.modules\.transformer",
    )

    import torch
    from src.models import AgnosticModelInferNoSmoother, AncestryLevelConvSmoother

    base_dir = args.model_dir / "Basemodel_Mamba_4096"
    smoother_dir = args.model_dir / "Smoother_512_Kernel_8192"
    base_args = load_pickle(base_dir / "args.pckl")
    base_model = AgnosticModelInferNoSmoother(base_args)
    base_state = torch.load(base_dir / "best_model.pth", map_location="cpu")
    base_model.load_state_dict(base_state)

    # The current upstream runner constructs the smoother from its published
    # architecture and uses args.pckl for inference window settings.
    load_pickle(smoother_dir / "args.pckl")
    smoother_model = AncestryLevelConvSmoother(
        in_planes=4, planes=32, kernel_size=8192
    )
    smoother_state = torch.load(smoother_dir / "best_model.pth", map_location="cpu")
    smoother_model.load_state_dict(smoother_state)
    print(
        "AS3 MODEL CHECK PASSED: "
        f"base_parameters={sum(p.numel() for p in base_model.parameters())} "
        f"smoother_parameters={sum(p.numel() for p in smoother_model.parameters())}"
    )


if __name__ == "__main__":
    main()
