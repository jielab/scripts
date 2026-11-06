"""Executed in the actual R/reticulate process after ems120.py is loaded."""
import tempfile
from pathlib import Path

import numpy as np
import torch
from sklearn.model_selection import cross_val_score
from xgboost import XGBRegressor
from transformers import AutoTokenizer, AutoModelForSequenceClassification, Trainer, TrainingArguments

# source_python publishes the analysis helper functions in __main__.
import __main__ as ems

training = ems._validate_training_df(ems.ems120_preflight_training, require_cv=True)
print(f"PASS: training labels; rows={len(training)}, classes={training[ems.LABEL_COL].nunique()}, stratified CV=5", flush=True)
del training, ems.ems120_preflight_training

assert ems.check_phone_gpu()
x = np.arange(48, dtype=np.float32).reshape(24, 2)
y = np.linspace(0, 10, 24, dtype=np.float32)
candidate = XGBRegressor(n_estimators=2, max_depth=2, n_jobs=ems.CPU_CORES, tree_method="hist", device="cuda", verbosity=0)
scores = cross_val_score(candidate, x, y, cv=3, scoring="neg_mean_squared_error", error_score="raise")
assert np.isfinite(scores).all()
print("PASS: XGBoost CUDA and scikit-learn cross-validation interoperability", flush=True)

model_root = Path(ems.BASE_MODEL_DIR)
assert (model_root / "config.json").is_file()
tokenizer = AutoTokenizer.from_pretrained(model_root, local_files_only=True)
model = AutoModelForSequenceClassification.from_pretrained(model_root, num_labels=2, local_files_only=True).to("cuda")
# Match the production training batch and maximum length with entirely
# synthetic Chinese strings. Only one optimizer step; never save the model.
encoding = tokenizer(["测试文本：用于环境验证。" * 16] * 16, truncation=True, padding="max_length", max_length=128)
dataset = ems.DxDataset(encoding, [i % 2 for i in range(16)])
with tempfile.TemporaryDirectory(prefix="ems120-hfl-check-") as out:
    args = TrainingArguments(output_dir=out, per_device_train_batch_size=16, max_steps=1,
                             save_strategy="no", logging_strategy="no", fp16=True,
                             report_to=[], disable_tqdm=True, seed=12010)
    result = Trainer(model=model, args=args, train_dataset=dataset).train()
    assert np.isfinite(result.training_loss)
    model.eval()
    inputs = {k: torch.tensor(v[:2], device="cuda") for k, v in encoding.items()}
    with torch.no_grad():
        logits = model(**inputs).logits
    assert tuple(logits.shape) == (2, 2) and torch.isfinite(logits).all()
print("PASS: offline HFL weights/tokenizer, FP16 Trainer step (batch 16, length 128) and GPU inference", flush=True)
del model, dataset, encoding, inputs, logits
torch.cuda.empty_cache()
