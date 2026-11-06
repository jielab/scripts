"""Independent trait jobs share immutable inputs and keep exclusive writers."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import threading

import pytest

from dosage_storage import dosage_lock, prsformer_cache_lock


def test_two_readers_exclude_input_replacement_and_release(tmp_path):
    config = dict(cache_root=str(tmp_path / "blocks"))
    cache = tmp_path / "genotype"
    with dosage_lock(config), dosage_lock(config), prsformer_cache_lock(cache), prsformer_cache_lock(cache):
        with pytest.raises(BlockingIOError):
            with dosage_lock(config, exclusive=True):
                pass
        with pytest.raises(BlockingIOError):
            with prsformer_cache_lock(cache, exclusive=True, wait=False):
                pass
    with dosage_lock(config, exclusive=True), prsformer_cache_lock(cache, exclusive=True, wait=False):
        pass


def test_prsformer_dispatch_allows_two_independent_trait_workers(tmp_path, monkeypatch):
    import prsformer_resume
    from abm_core import cold_start
    barrier = threading.Barrier(2)
    config = dict(cache_root=str(tmp_path / "blocks"))
    monkeypatch.setattr(cold_start, "benchmark_config", lambda args: {})
    monkeypatch.setattr(cold_start, "benchmark_device", lambda device: "cuda")
    monkeypatch.setattr(cold_start, "ensure_inputs", lambda *args: None)

    def worker(code_root, options, configuration, device, trait):
        # The tested dispatcher still owns real kinship and trait output locks.
        # Only GPU work is replaced; independent jobs must enter together.
        with dosage_lock(config), prsformer_cache_lock(tmp_path / "genotype"):
            barrier.wait(timeout=5)
            return trait

    monkeypatch.setattr(prsformer_resume, "_run_trait", worker)
    def dispatch(trait):
        prsformer_resume.dispatch(Path(__file__).parents[1],
            ["--analysis-root", str(tmp_path), "--traits", trait, "--device", "cuda"])
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(dispatch, ["height", "ldl"]))
    assert (tmp_path / "prsformer/height/.resume.lock").is_file()
    assert (tmp_path / "prsformer/ldl/.resume.lock").is_file()
    assert not (tmp_path / "prsformer/.resume.lock").exists()


def test_two_trait_reconstruction_probes_can_read_shared_dosage(tmp_path):
    import block_inputs as b
    from abm_core.grid_inputs import locate_blocks
    from test_block_producer import fixture
    config_path, config = fixture(tmp_path, "ALL.37.bed")
    shutil.copytree(tmp_path / "weights/height", tmp_path / "weights/ldl")
    shutil.copytree(tmp_path / "scores/height", tmp_path / "scores/ldl")
    for trait in ("height", "ldl"):
        b.prepare_trait(config_path, trait, device="cpu")
    outer = dict(block_config=str(config_path), execution={"device": "cpu"})
    barrier = threading.Barrier(2)
    def read(trait):
        barrier.wait(timeout=5)
        array, _, _, receipt, _ = locate_blocks(outer, trait, allow_build=False)
        array._mmap.close()
        return receipt["trait"]
    with ThreadPoolExecutor(max_workers=2) as pool:
        assert set(pool.map(read, ["height", "ldl"])) == {"height", "ldl"}
