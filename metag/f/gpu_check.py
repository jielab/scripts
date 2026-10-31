#!/usr/bin/env python3
from __future__ import annotations
import re, shutil, subprocess
from common import DIR_OUT, PART_COMPRESSED_GIB, atomic_json

def capture(cmd):
    return subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT).strip()

def report(result):
    print(f"GPU detection: {'CUDA available' if result['cuda_gpu'] else 'no usable NVIDIA GPU; CPU fallback'}")
    for i, gpu in enumerate(result.get("gpus", [])):
        print(f"  GPU {i}: {gpu['name']} (compute capability {gpu['compute_capability']})")
        print(f"    VRAM: {gpu['memory_total_mib']/1024:.1f} GiB total, "
              f"{gpu['memory_used_mib']/1024:.1f} GiB used, "
              f"{gpu['memory_free_mib']/1024:.1f} GiB free")
    print(f"CUDA Toolkit: {result.get('cuda_toolkit', 'not found')}")
    if "recommended_cuda_arch" in result:
        print(f"CUDA build target: {result['recommended_cuda_arch']}")
        print("Toolkit supports Blackwell: " + ("yes" if result.get("cuda_toolkit_ok_for_blackwell") else "no"))
    for key in ("nvidia_smi_error", "nvcc_error"):
        if key in result:
            print(f"  {key}: {result[key]}")
    if "recommended_partition_gib_now" in result:
        print(f"Estimated VRAM budget per partition now: {result['recommended_partition_gib_now']} GiB "
              "(not compressed FASTA size)")
    print(f"Configured partition input: {PART_COMPRESSED_GIB:g} GiB compressed FASTA "
          "(overridable with --part-gib)")


def main():
    result = {"cuda_gpu": False, "gpu_backend": "cpu"}
    if shutil.which("nvidia-smi"):
        q = ["nvidia-smi", "--query-gpu=name,memory.total,memory.used,memory.free,compute_cap",
             "--format=csv,noheader,nounits"]
        try:
            rows = []
            for line in capture(q).splitlines():
                name,total,used,free,cc = [x.strip() for x in line.split(",",4)]
                rows.append(dict(name=name, memory_total_mib=int(total),
                                 memory_used_mib=int(used), memory_free_mib=int(free),
                                 compute_capability=cc))
            if rows:
                result.update(cuda_gpu=True, gpu_backend="cuda", gpus=rows)
        except Exception as e:
            result["nvidia_smi_error"] = str(e)

    if shutil.which("nvcc"):
        try:
            out = capture(["nvcc","--version"])
            result["nvcc"] = out
            m = re.search(r"release\s+(\d+)\.(\d+)", out)
            if m:
                result["cuda_toolkit"] = f"{m.group(1)}.{m.group(2)}"
                result["cuda_toolkit_ok_for_blackwell"] = (int(m.group(1)),int(m.group(2))) >= (12,8)
        except Exception as e:
            result["nvcc_error"] = str(e)

    if result.get("gpus"):
        minfree = min(x["memory_free_mib"] for x in result["gpus"])
        result["recommended_partition_gib_now"] = max(0, min(12, int(minfree * 0.65 / 1024)))
        # RTX 50 series / CC 12.0
        if any(x.get("compute_capability") == "12.0" for x in result["gpus"]):
            result["recommended_cuda_arch"] = "sm_120"

    out = DIR_OUT / "gpu" / "gpu.json"
    atomic_json(out, result)
    report(result)
    print(f"GPU report saved: {out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
