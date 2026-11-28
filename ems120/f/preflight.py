"""Synthetic checks in the real allocated R/reticulate process; no production fit."""
import os
import tempfile
from pathlib import Path

import numpy as np
import torch
from sklearn.model_selection import cross_val_score
from transformers import AutoTokenizer, AutoModelForSequenceClassification, Trainer, TrainingArguments, set_seed

# source_python publishes the analysis helpers in __main__.
import __main__ as ems

phase = os.environ['EMS120_PHASE']
device = ems.DEVICE
training = ems._validate_training_df(ems.ems120_preflight_training, require_cv=True)
print(f"PASS: training labels; rows={len(training)}, classes={training[ems.LABEL_COL].nunique()}, stratified CV=5", flush=True)
del training, ems.ems120_preflight_training

if phase == 'phone':
    # Check all configured comparators without changing production candidates.
    candidates = ems._candidate_models()
    if ems.PHONE_RUN == 'gpu':
        assert ems.check_phone_gpu() and 'xgboost' in candidates
    x = np.arange(60, dtype=np.float32).reshape(30, 2)
    y = np.linspace(0, 10, 30, dtype=np.float32)
    for name, candidate in candidates.items():
        params = candidate.get_params()
        if 'n_estimators' in params:
            candidate.set_params(n_estimators=2)
        scores = cross_val_score(candidate, x, y, cv=3, scoring='neg_mean_squared_error', error_score='raise')
        assert np.isfinite(scores).all(), name
    print(f"PASS: phone comparator/scikit-learn synthetic CV; XGBoost processor={ems.PHONE_RUN}", flush=True)
else:
    assert phase in {'geo', 'dx'}
    model_root = Path(ems.BASE_MODEL_DIR)
    assert (model_root / 'config.json').is_file()
    set_seed(ems.DX_CV_SEED)
    tokenizer = AutoTokenizer.from_pretrained(model_root, local_files_only=True)
    model = AutoModelForSequenceClassification.from_pretrained(model_root, num_labels=2, local_files_only=True).to(device)
    length = 96 if phase == 'geo' else 128
    encoding = tokenizer(['测试文本：用于环境验证。' * 16] * 16, truncation=True, padding='max_length', max_length=length)
    dataset = ems.DxDataset(encoding, [i % 2 for i in range(16)])
    with tempfile.TemporaryDirectory(prefix=f'ems120-{phase}-check-') as out:
        args = TrainingArguments(output_dir=out, per_device_train_batch_size=16, max_steps=1,
                                 save_strategy='no', logging_strategy='no', use_cpu=(device == 'cpu'),
                                 fp16=(device == 'cuda'), report_to=[], disable_tqdm=True,
                                 seed=ems.DX_CV_SEED, data_seed=ems.DX_CV_SEED)
        result = Trainer(model=model, args=args, train_dataset=dataset).train()
        assert np.isfinite(result.training_loss)
        model.eval()
        inputs = {k: torch.tensor(v[:2], device=device) for k, v in encoding.items()}
        with torch.no_grad():
            logits = model(**inputs).logits
        assert tuple(logits.shape) == (2, 2) and torch.isfinite(logits).all()
    print(f'PASS: offline HFL {phase} one-step training/inference, device={device}, fp16={device == "cuda"}', flush=True)
    del model, dataset, encoding, inputs, logits
    if device == 'cuda':
        torch.cuda.empty_cache()
