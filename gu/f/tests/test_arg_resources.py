"""Run under WSL: python3 f/tests/test_arg_resources.py [--live]. No genetic data."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

GU = Path(__file__).resolve().parents[2]


def run(script, env=None):
    return subprocess.run(["bash", "-c", script], text=True, capture_output=True,
                          env=env, timeout=30)


with tempfile.TemporaryDirectory(prefix="gu-resource-test-") as temp:
    root = Path(temp)
    (root / "f").mkdir()
    mock = root / "bin"
    mock.mkdir()
    for name, body in {
        "systemctl": '#!/bin/bash\necho "${TEST_SYSTEMD_STATE:-running}"\n',
        "systemd-run": '#!/bin/bash\nprintf "%s\\n" "$@"\nexit "${TEST_SCOPE_STATUS:-0}"\n',
    }.items():
        path = mock / name
        path.write_text(body)
        path.chmod(0o755)
    env = dict(os.environ, PATH=f"{mock}:{os.environ['PATH']}",
               ARG_MEMORY_LIMIT_GB="32", ARG_MEMORY_SWAP_GB="4")
    # Exercise the public entry point, including argument consumption.
    def public(*args, overrides=None):
        return subprocess.run(["bash", str(GU / "arg.sh"), *args],
                              env=dict(env, **(overrides or {})), text=True,
                              capture_output=True, timeout=30)

    result = public("build", "--dir-gen", "/unused path", "--method", "needle",
                    "--format", "trace", "--threads", "8", "--jobs", "4")
    assert result.returncode == 0, result
    assert "MemoryMax=32G" in result.stdout and "MemorySwapMax=4G" in result.stdout
    assert "--jobs\n1\n" in result.stdout and "/unused path\n" in result.stdout
    assert "--threads\n8\n" in result.stdout and "reducing build" in result.stderr
    result = public("prep_gen", "--jobs", "4", "--memory-limit-gb", "24",
                    "--memory-swap-gb", "0")
    assert result.returncode == 0 and "--jobs\n4\n" in result.stdout
    assert "concurrency is limited to 1" not in result.stderr
    assert "MemoryMax=24G" in result.stdout and "MemorySwapMax=0G" in result.stdout
    assert "--memory-limit-gb" not in result.stdout
    result = public("check")
    assert "--action\ncheck\n" in result.stdout
    for args in [("--memory-limit-gb", "0"), ("--memory-limit-gb", "bad"),
                 ("--memory-swap-gb", "-1"), ("--jobs", "0"), ("--jobs",)]:
        result = public("build", *args)
        assert result.returncode == 2 and not result.stdout, result
    result = public("build", overrides={"TEST_SYSTEMD_STATE": "offline"})
    assert result.returncode == 2 and "refusing" in result.stderr and not result.stdout
    result = public("build", overrides={"TEST_SCOPE_STATUS": "137"})
    assert result.returncode == 137 and "does not identify the cause" in result.stderr
    print("PASS: public CLI caps build concurrency, preserves arguments, fails closed and propagates errors")

    if "--live" in sys.argv:
        driver = root / "f" / "arg.sh"
        driver.write_text('''#!/bin/bash
set -eu
group=$(cut -d: -f3 /proc/self/cgroup)
cat "/sys/fs/cgroup${group}/memory.max"
cat "/sys/fs/cgroup${group}/memory.swap.max"
cat "/sys/fs/cgroup${group}/memory.oom.group"
bash -c 'cat /proc/self/cgroup'
''')
        prefix = f'ROOT={str(root)!r}; source {str(GU / "f/arg_resources.sh")!r}; '
        result = run(prefix + 'arg_run_guarded build --jobs 4')
        assert result.returncode == 0, result
        assert "34359738368\n4294967296\n1\n" in result.stdout, result
        assert ".scope" in result.stdout
        print("PASS: live scope enforces 32 GiB RAM / 4 GiB swap and children inherit it")

        # A deliberately tiny test-only cap verifies OOM containment, without
        # allocating a large workload or touching real pipeline inputs.
        driver.write_text("#!/bin/bash\npython3 -c 'x=bytearray(128*1024*1024)'\n")
        substitute = '''systemd-run() {
  local -a values=(); local a
  for a in "$@"; do
    case "$a" in MemoryMax=*) a=MemoryMax=32M;; MemorySwapMax=*) a=MemorySwapMax=0;; esac
    values+=("$a")
  done
  command systemd-run "${values[@]}"
}
'''
        result = run(prefix + substitute + '''
if arg_run_guarded build; then exit 99; else echo "PARENT_SURVIVED:$?"; fi
''')
        assert result.returncode == 0 and "PARENT_SURVIVED:" in result.stdout, result
        assert "ARG exited with status" in result.stderr, result
        print("PASS: bounded OOM kills the test workload; launcher outside the scope survives")
