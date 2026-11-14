"""Conservative single-node CPU budgets from the user's QiMing v2.0 inventory.

These are planning budgets (GiB), not live availability or site policy. Keep
headroom for the OS. Unknown explicit queues remain available for other sites.
No queue discovery, submission, numerical imports or data access occurs here.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, InvalidOperation, ROUND_CEILING
import json
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
    "9654": (192, 704, 720),
}

# Physical RAM, not the smaller application budgets above. Admin's 2026-11-05
# rule reserves quarters of a CPU node based on max(CPU fraction, RAM fraction).
NODE_MEMORY_GIB = {"33": 192, "38": 512, "75m": 192, "9654": 768,
                   "v3-64": 64, "73x": 512, "63": 512, "52": 512,
                   "ot38": 2560, "2t50c": 2048, "1t75c": 1024,
                   "1t81c": 1024, "1t88c": 1024, "v4-3t": 3072,
                   "v3-6t": 6144, "sk-6t-104": 6144}
# These two queues were explicitly named by Admin as allowing exact CPU counts.
EXACT_CORE_QUEUES = {"v4-3t", "v3-6t"}
POLICY_VERSION = "qiming-admin-20261105-quarter-node"


def cpu_slots(queue: str, cores: int, memory_gib: float, hours: float) -> tuple[str, int]:
    """Return queue and allocated slots for a NEW single-node CPU submission.

    Scientific thread counts remain independent. GPU/unknown queues are not
    assigned an invented node shape. Call only while planning/submitting HPC,
    never when loading an existing frozen plan or executing a local workflow.
    """
    try:
        c, m, h = (Decimal(str(x)) for x in (cores, memory_gib, hours))
    except InvalidOperation as exc:
        raise ValueError("CPU cores, memory and hours must be numbers") from exc
    if any(not x.is_finite() or x <= 0 for x in (c, m, h)) or c != c.to_integral_value():
        raise ValueError("CPU cores, memory and hours must be positive; cores must be integral")
    queue = cpu_queue(queue, int(c), float(m), float(h))
    if queue not in NODE_MEMORY_GIB or queue in EXACT_CORE_QUEUES:
        return queue, int(c)
    node_cores = CPU_LIMITS[queue][0]
    fraction = max(c / node_cores, m / NODE_MEMORY_GIB[queue])
    quarters = int((fraction * 4).to_integral_value(rounding=ROUND_CEILING))
    allocated = int((Decimal(node_cores) * quarters / 4).to_integral_value(rounding=ROUND_CEILING))
    return queue, allocated


def main():
    p = argparse.ArgumentParser(description="QiMing Admin CPU/RAM quarter-node recommendation (2026-11-05)")
    p.add_argument('--queue', default='33')
    p.add_argument('--cores', type=int, required=True)
    p.add_argument('--memory-gib', type=float, required=True)
    duration = p.add_mutually_exclusive_group()
    duration.add_argument('--hours', type=float, default=24)
    duration.add_argument('--walltime', help='LSF HHH:MM walltime, including minutes')
    p.add_argument('--format', choices=['json', 'shell', 'bsub'], default='json')
    a = p.parse_args()
    try:
        hours = a.hours
        if a.walltime is not None:
            parts = a.walltime.split(':')
            if len(parts) != 2 or not all(x.isdigit() for x in parts) or not 0 <= int(parts[1]) < 60:
                raise ValueError('walltime must be HHH:MM, with minutes in 00..59')
            hours = Decimal(parts[0]) + Decimal(parts[1]) / 60
        q, n = cpu_slots(a.queue, a.cores, a.memory_gib, hours)
    except ValueError as exc:
        p.error(str(exc))
    if a.format == 'shell':
        print(q, n)
    elif a.format == 'bsub':
        print(f'-q {q} -n {n} -R "span[hosts=1] rusage[mem={a.memory_gib:g}GB/host]"')
    else:
        print(json.dumps(dict(queue=q, requested_cores=a.cores, allocated_cores=n,
                              memory_gib=a.memory_gib, policy=POLICY_VERSION)))


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


if __name__ == '__main__':
    main()
