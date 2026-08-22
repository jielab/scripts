#!/usr/bin/env python3
"""Download/cache PubMedBERT/BiomedBERT model files before training."""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import annotations

import argparse
from pathlib import Path

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model_name", default="microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext")
    ap.add_argument("--fallback_model_name", default="microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext")
    ap.add_argument("--cache_dir", default=None)
    args = ap.parse_args()

    try:
        from huggingface_hub import snapshot_download
    except Exception as e:
        raise RuntimeError("Please install huggingface_hub first: pip install huggingface_hub") from e

    for name in [args.model_name, args.fallback_model_name]:
        try:
            path = snapshot_download(repo_id=name, cache_dir=args.cache_dir)
            print(f"Cached {name} at {path}")
            return
        except Exception as e:
            print(f"WARNING: failed to download/cache {name}: {e}")
    raise RuntimeError("Could not download either PubMedBERT model. Check network/Hugging Face access.")

if __name__ == "__main__":
    main()
