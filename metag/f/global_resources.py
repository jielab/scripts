"""Conservative Linux/WSL memory preflight; never change host memory settings.

The database budget limits Kraken's hash table, not its process peak RSS.
These checks are admission checks; concurrent processes and build buffers can
still consume memory after the check. We deliberately do not auto-shrink a DB.
"""
from __future__ import annotations

import math
from pathlib import Path

GIB = 1 << 30


def _integer(path):
    try:
        value = Path(path).read_text().strip()
    except FileNotFoundError:
        return None
    except OSError as exc:
        raise RuntimeError(f"Cannot read cgroup memory value: {path}") from exc
    if value == "max":
        return None
    try:
        number = int(value)
    except ValueError as exc:
        raise RuntimeError(f"Invalid cgroup memory value at {path}: {value!r}") from exc
    if number < 0:
        raise RuntimeError(f"Invalid negative cgroup memory value at {path}: {number}")
    return number


def _cgroup_paths(proc_root, cgroup_root):
    """Include the exposed root plus the process group and its ancestors.

    Container cgroup namespaces may report '/' or a host-relative membership
    which is not visible below the mounted root. The root check covers both.
    """
    roots = {cgroup_root, cgroup_root / "memory"}
    memberships = []
    try:
        for line in (proc_root / "self/cgroup").read_text().splitlines():
            _hierarchy, controllers, relative = line.split(":", 2)
            if not controllers or "memory" in controllers.split(","):
                memberships.append((controllers, relative))
    except (OSError, ValueError):
        pass
    candidates = set(roots)
    for controllers, relative in memberships:
        # '..' is possible in a cgroup namespace; never escape an exposed root.
        parts = Path(relative.lstrip("/")).parts
        if ".." in parts:
            continue
        for root in (cgroup_root, cgroup_root / "memory"):
            current = root.joinpath(*parts)
            while current == root or root in current.parents:
                candidates.add(current)
                if current == root:
                    break
                current = current.parent
    return sorted(candidates, key=str)


def memory_status(*, proc_root=Path("/proc"), cgroup_root=Path("/sys/fs/cgroup")):
    proc_root, cgroup_root = Path(proc_root), Path(cgroup_root)
    values = {}
    try:
        for line in (proc_root / "meminfo").read_text().splitlines():
            key, raw = line.split(":", 1)
            fields = raw.split()
            values[key] = int(fields[0]) * (1024 if len(fields) > 1 and fields[1] == "kB" else 1)
    except (OSError, ValueError, IndexError) as exc:
        raise RuntimeError("Cannot read Linux /proc/meminfo; memory preflight is required") from exc
    total = values.get("MemTotal", 0)
    available = values.get("MemAvailable")
    if available is None:
        available = (values.get("MemFree", 0) + values.get("Buffers", 0) +
                     values.get("Cached", 0) + values.get("SReclaimable", 0) -
                     values.get("Shmem", 0))
    if total <= 0 or available < 0:
        raise RuntimeError("Invalid total/available memory in /proc/meminfo")
    effective = min(total, available)
    limits = []
    for folder in _cgroup_paths(proc_root, cgroup_root):
        for version, limit_name, usage_name in (
            (2, "memory.max", "memory.current"),
            (1, "memory.limit_in_bytes", "memory.usage_in_bytes"),
        ):
            limit = _integer(folder / limit_name)
            # Linux v1 represents an unlimited group with a huge sentinel.
            if limit is None or limit >= (1 << 60):
                continue
            usage = _integer(folder / usage_name)
            if usage is None:
                raise RuntimeError(f"Cannot read usage for finite cgroup memory limit: {folder}")
            remaining = max(0, limit - usage)
            effective = min(effective, remaining)
            limits.append({"version": version, "path": str(folder),
                           "limit_bytes": limit, "usage_bytes": usage,
                           "remaining_bytes": remaining})
    return {"mem_total_bytes": total, "mem_available_bytes": available,
            "effective_available_bytes": effective, "cgroups": limits,
            "note": "Preflight only; hash budget is not a peak process memory limit."}


def require_memory(required_bytes, reserve_gib=16, context="global Kraken 2"):
    if (not isinstance(required_bytes, (int, float)) or not math.isfinite(required_bytes)
            or required_bytes < 0 or not math.isfinite(reserve_gib) or reserve_gib < 0):
        raise ValueError("Memory requirement and reserve must be finite and nonnegative")
    required = math.ceil(required_bytes)
    reserve = math.ceil(reserve_gib * GIB)
    state = memory_status()
    state.update(required_bytes=required, reserve_bytes=reserve, context=context)
    if state["effective_available_bytes"] < required + reserve:
        raise RuntimeError(
            f"Insufficient Linux/WSL/cgroup available memory for {context}: "
            f"{state['effective_available_bytes']/GIB:.1f} GiB available, "
            f"{required/GIB:.1f} GiB database budget/files + {reserve/GIB:.1f} GiB reserve required. "
            "Close competing tasks or deliberately choose a smaller --db-gib when building. "
            "A smaller database changes minimizer sampling and may change sensitivity; "
            "no automatic budget reduction is performed."
        )
    return state
