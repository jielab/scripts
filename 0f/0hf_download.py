#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#%%
import os
from huggingface_hub import snapshot_download

# Configure a mirror for WSL or restricted network environments.
os.environ["HF_ENDPOINT"] = "https://hf-mirror.com"
local_dir = "/mnt/d/data.BIG/AI_model_files/hfl"

print("正在开始下载...")
snapshot_download(
    token="hf_XXX", 
    repo_id="hfl/chinese-macbert-large", 
    repo_type="model", 
    local_dir=local_dir,
    local_dir_use_symlinks=False
)
print(f"下载完成！文件已存至: {local_dir}")