# 🚩 Imports, paths, and shared inputs
# %%

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import os
from huggingface_hub import snapshot_download

# Configure a mirror for WSL or restricted network environments.
os.environ["HF_ENDPOINT"] = "https://hf-mirror.com"
local_dir = configured_path('DIR_BIG', 'AI_model/hfl')

print("正在开始下载...")
snapshot_download(
	token="hf_XXX",
	repo_id="hfl/chinese-macbert-large",
	repo_type="model",
	local_dir=local_dir,
	local_dir_use_symlinks=False,
)
print(f"下载完成！文件已存至: {local_dir}")
