#!/usr/bin/env python3
"""Validate and report the isolated ArchaicSeeker3 Python/Torch/CUDA stack."""
from __future__ import annotations

import importlib
import os
import re
import subprocess
import sys
import tempfile


def major_minor(value: str | None) -> str:
    match = re.search(r"(\d+)\.(\d+)", value or "")
    return f"{match.group(1)}.{match.group(2)}" if match else ""


def gpu_compute_capability() -> str:
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader,nounits"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
        return major_minor(result.stdout.splitlines()[0]) if result.stdout.strip() else ""
    except Exception:
        return ""


def main() -> None:
    expected_python = os.environ.get("AS3_EXPECT_PYTHON") or os.environ.get("AS3_PYTHON_VERSION") or "3.9"
    compute_cap = os.environ.get("AS3_GPU_COMPUTE_CAP", "") or gpu_compute_capability()
    try:
        blackwell = int((compute_cap or "0").split(".", 1)[0]) >= 12
    except ValueError:
        blackwell = False
    profile = "blackwell" if blackwell else "upstream"
    expected_torch = os.environ.get("AS3_EXPECT_TORCH") or os.environ.get("AS3_TORCH_VERSION") or ("2.8.0" if blackwell else "2.4.0")
    expected_cuda = os.environ.get("AS3_EXPECT_CUDA") or os.environ.get("AS3_CUDA_VERSION") or ("12.8" if blackwell else "12.1")
    require_cuda = os.environ.get("AS3_REQUIRE_CUDA", "1") == "1"

    errors: list[str] = []
    actual_python = f"{sys.version_info.major}.{sys.version_info.minor}"
    if actual_python != expected_python:
        errors.append(f"Python {actual_python} != expected {expected_python}")

    torch = None
    try:
        torch = importlib.import_module("torch")
    except Exception as exc:
        errors.append(f"torch import failed: {exc}")

    # Import the full inference dependency chain here.  Checking only Torch and
    # mamba-ssm allowed an environment with an incompatible NumPy/Dask pair (or
    # without torchvision) to pass this check and then fail every scheduled
    # chromosome task after GPU allocation.
    for module in (
        "mamba_ssm",
        "causal_conv1d",
        "numpy",
        "pandas",
        "pysam",
        "h5py",
        "dask.array",
        "allel",
        "sklearn",
        "torchvision",
        "torchmetrics",
        "transformers",
        "einops",
        "cv2",
    ):
        try:
            importlib.import_module(module)
        except Exception as exc:
            detail = str(exc).splitlines()[0] if str(exc) else type(exc).__name__
            errors.append(f"{module} import failed: {detail}")

    h5py = sys.modules.get("h5py")
    if h5py is not None:
        built_hdf5 = tuple(getattr(h5py.version, "hdf5_built_version_tuple", ()))
        runtime_hdf5 = tuple(getattr(h5py.version, "hdf5_version_tuple", ()))
        if built_hdf5[:2] != runtime_hdf5[:2]:
            errors.append(
                "h5py HDF5 runtime "
                f"{h5py.version.hdf5_version} != built "
                f"{'.'.join(map(str, built_hdf5))}; the installed HDF5 runtime is inconsistent"
            )
        h5_path = ""
        try:
            with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as handle:
                h5_path = handle.name
            with h5py.File(h5_path, "w") as handle:
                handle.create_dataset("gu_health", data=[1, 2, 3])
            with h5py.File(h5_path, "r") as handle:
                if list(handle["gu_health"][:]) != [1, 2, 3]:
                    raise RuntimeError("round-trip data mismatch")
        except Exception as exc:
            errors.append(f"h5py HDF5 read/write round-trip failed: {exc}")
        finally:
            if h5_path:
                try:
                    os.unlink(h5_path)
                except FileNotFoundError:
                    pass

    actual_torch = "missing"
    actual_cuda = "missing"
    cuda_available = False
    gpu_name = "none"
    torch_cap = ""
    if torch is not None:
        actual_torch = str(torch.__version__).split("+", 1)[0]
        actual_cuda = str(torch.version.cuda or "")
        cuda_available = bool(torch.cuda.is_available())
        if cuda_available:
            gpu_name = torch.cuda.get_device_name(0)
            torch_cap = ".".join(map(str, torch.cuda.get_device_capability(0)))
        if actual_torch != expected_torch:
            errors.append(f"torch {actual_torch} != expected {expected_torch}")
        if major_minor(actual_cuda) != major_minor(expected_cuda):
            errors.append(f"torch CUDA {actual_cuda or 'none'} != expected {expected_cuda}")
        if require_cuda and not cuda_available:
            errors.append("torch.cuda.is_available() is false")

    print(
        "AS3 HEALTH "
        f"profile={profile} python={actual_python} torch={actual_torch} "
        f"torch_cuda={actual_cuda or 'none'} cuda_available={str(cuda_available).lower()} "
        f"gpu={gpu_name} compute_cap={torch_cap or compute_cap or 'unknown'} "
        f"expected_python={expected_python} expected_torch={expected_torch} expected_cuda={expected_cuda}"
    )
    if errors:
        for error in errors:
            print(f"AS3 HEALTH ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
if __name__ == "__main__":
    main()
