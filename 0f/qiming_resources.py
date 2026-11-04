"""Conservative single-node CPU budgets from the user's QiMing v2.0 inventory.

These are planning budgets (GiB), not live availability or site policy. Keep
headroom for the OS. Unknown explicit queues remain available for other sites.
No queue discovery, submission, numerical imports or data access occurs here.
"""
from __future__ import annotations
import os

# queue: (cores per node, recommended maximum job RAM GiB, walltime hours)
CPU_LIMITS = {
    "33": (40, 160, 720), "38": (64, 448, 720),
    "ot38": (64, 1792, 336), "2t50c": (64, 1792, 336),
    "1t75c": (64, 896, 336), "1t88c": (96, 896, 336),
    "1t81c": (120, 896, 336), "73x": (128, 448, 336),
    "52": (48, 448, 336), "63": (128, 448, 336),
    "v4-3t": (72, 2560, 336), "sk-6t-104": (104, 5120, 336),
    "v3-6t": (144, 5120, 336), "75m": (48, 160, 720),
    "v3-64": (24, 56, 1440),
}


def cpu_queue(queue: str, cores: int, memory_gib: float, hours: float) -> str:
    """Resolve auto, or reject a known CPU queue with an oversized request."""
    def fits(name):
        c, m, h = CPU_LIMITS[name]
        return cores <= c and memory_gib <= m and hours <= h

    if queue == "auto":
        candidates = [os.getenv("HPC_CPU_QUEUE", "33"),
                      os.getenv("HPC_MEDIUM_QUEUE", "38"),
                      os.getenv("HPC_HIGHMEM_QUEUE", "1t81c"),
                      "73x", "2t50c", "v4-3t",
                      os.getenv("HPC_HUGEMEM_QUEUE", "v3-6t")]
        for candidate in dict.fromkeys(candidates):
            if candidate in CPU_LIMITS and fits(candidate):
                return candidate
        raise ValueError("No configured QiMing single-node CPU queue fits "
                         f"cores={cores}, memory={memory_gib:g} GiB, hours={hours:g}; "
                         "reduce within-job concurrency or select a verified explicit queue")
    if queue in CPU_LIMITS and not fits(queue):
        c, m, h = CPU_LIMITS[queue]
        raise ValueError(f"Queue {queue}: request {cores} cores/{memory_gib:g} GiB/{hours:g} h "
                         f"exceeds the conservative single-node budget {c} cores/{m} GiB/{h} h. "
                         "Use --queue auto, choose a larger queue, or reduce --jobs/--threads.")
    return queue
