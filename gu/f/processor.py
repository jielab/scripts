"""Canonical, method-specific processor selection for GU planners and workers."""
from __future__ import annotations
import os

METHODS = ("ibdmix", "trace", "phyml")


def resolve(value: str | None) -> str:
    value = "cpu" if value is None else value
    if value not in {"cpu", "gpu", "auto"}:
        raise ValueError("processor must be cpu, gpu or auto")
    return "cpu" if value == "auto" else value


def selection(method: str, values: dict[str, str | None]) -> dict[str, str]:
    eligible = {"trace"} if method == "arg-trace" else {method}
    for name, value in values.items():
        if name not in METHODS:
            raise ValueError(f"unknown processor method: {name}")
        if value is not None and name not in eligible:
            raise ValueError(f"--{name}-run is not valid for {method}")
    return {name: resolve(values.get(name)) for name in METHODS}


def from_argv(method: str, argv: list[str]) -> dict[str, str]:
    values = {}
    for i, arg in enumerate(argv):
        if arg == "--GPU" or arg.startswith("--GPU=") or arg in {"trace-GPU", "phyml-GPU"}:
            raise ValueError("Use the method name and --ibdmix-run, --trace-run or --phyml-run cpu|gpu|auto")
        for name in METHODS:
            flag = f"--{name}-run"
            if arg == flag or arg.startswith(flag + "="):
                if name in values:
                    raise ValueError(f"Duplicate option: {flag}")
                if arg == flag:
                    if i + 1 >= len(argv) or argv[i + 1].startswith("--"):
                        raise ValueError(f"{flag} requires cpu, gpu or auto")
                    values[name] = argv[i + 1]
                else:
                    values[name] = arg.split("=", 1)[1]
    return selection(method, values)


def uses_gpu(plan: dict, stage: str) -> bool:
    method = plan["method"]
    if stage == "trace" and method == "arg-trace":
        return plan["processors"]["trace"] == "gpu"
    return stage == "worker" and method in METHODS and plan["processors"][method] == "gpu"


def validate_queue(queue: str, gpu: bool) -> None:
    """Reject known CPU/GPU conflicts; leave unknown site queues available."""
    from qiming_resources import CPU_LIMITS
    if gpu and queue in CPU_LIMITS:
        raise ValueError(f"Queue {queue} is in the configured CPU inventory; select a GPU queue for a GPU method")
    configured_gpu = {os.environ.get(name) for name in ("GU_GPU_QUEUE", "HPC_GPU_QUEUE")}
    if not gpu and queue not in CPU_LIMITS and queue in configured_gpu:
        raise ValueError(f"Queue {queue} is configured as a GPU queue; select a CPU queue for this CPU stage")
