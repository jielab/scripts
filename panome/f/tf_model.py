"""TabICLv2 fine-tuning with IPCW loss and fold-disjoint context labels.

The official sklearn fit() only prepares in-context learning. This adapter runs
actual AdamW gradient steps using the same raw Transformer forward as the
official fine-tuner, adding exact query loss weights and family exclusion.
The private numerical adapter/model API is pinned to tabicl==2.2.0.
"""
from pathlib import Path
import time
import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from common import dump, log
from survival import loss


def select_features(x, y, w, maximum):
    """Weighted univariate association, computed only on the build sample."""
    w = np.asarray(w, float)
    w = w/w.sum()
    mean_y = w @ y
    mean_x = w @ x
    cov = (w*(y-mean_y)) @ x
    variance = np.maximum(w @ np.square(x.astype(float))-mean_x**2, 1e-12)
    score = np.abs(cov)/np.sqrt(variance*max(mean_y*(1-mean_y), 1e-12))
    n = x.shape[1] if maximum == 0 else min(maximum, x.shape[1])
    selected = np.sort(np.argsort(-score, kind="stable")[:n])
    return selected, score


def context_indices(pool, y, w, size, rng):
    """Unique labeled context; preserve weighted class proportions, not 50:50."""
    pool = np.asarray(pool, int)
    size = min(size, len(pool))
    cases, controls = pool[y[pool] == 1], pool[y[pool] == 0]
    if min(len(cases), len(controls)) < 2 or size < 4:
        raise ValueError("Each context pool needs at least two known cases and controls")
    nc = int(round(size*w[cases].sum()/w[pool].sum()))
    nc = min(len(cases), size-2, max(2, size-len(controls), nc))
    choices = []
    for part, count in [(cases, nc), (controls, size-nc)]:
        choices.extend(rng.choice(part, count, replace=False, p=w[part]/w[part].sum()))
    return rng.permutation(choices)


def fold_labels(y, groups, count, seed):
    cv = (StratifiedKFold(count, shuffle=True, random_state=seed) if groups is None else
          StratifiedGroupKFold(count, shuffle=True, random_state=seed))
    folds = np.empty(len(y), int)
    for fold, (_, ix) in enumerate(cv.split(np.zeros(len(y)), y, groups)):
        folds[ix] = fold
    if groups is not None:
        frame = pd.DataFrame({"group": groups, "fold": folds})
        if frame.groupby("group").fold.nunique().max() > 1:
            raise ValueError("Family leakage across fine-tuning context folds")
    return folds


def _predictor(model, config, cx, cy, a):
    from tabicl import TabICLClassifier
    estimator = TabICLClassifier(n_estimators=1, norm_methods=["none"],
        feat_shuffle_method="none", class_shuffle_method="none", softmax_temperature=1.,
        device=a.device, use_amp=a.amp and a.device == "cuda", use_fa3=False,
        n_jobs=a.cores, kv_cache=True, allow_auto_download=False, random_state=a.seed)
    estimator.model_, estimator.model_config_, estimator.model_path_ = model, config, a.model_path
    # This is the same injection used by the official fine-tuning wrapper.
    estimator._load_model = lambda: None
    model.eval()
    return estimator.fit(cx, cy)


def predict(estimator, x, batch_size):
    rows = []
    last = time.monotonic()
    for start in range(0, len(x), batch_size):
        pred = estimator.predict_proba(x[start:start+batch_size])[:, 1]
        if not np.isfinite(pred).all():
            raise FloatingPointError("Nonfinite Transformer prediction")
        rows.append(pred)
        if time.monotonic()-last > 45:
            log("RUNNING", "TF predict", f"rows={start+len(pred)}/{len(x)}")
            last = time.monotonic()
    return np.concatenate(rows) if rows else np.empty(0)


def train_transformer(x, y, w, ids, groups, xv, yv, wv, a, out):
    from tabicl._model import TabICL
    from tabicl._sklearn.preprocessing import EnsembleGenerator
    torch.manual_seed(a.seed)
    if a.device == "cuda":
        torch.cuda.manual_seed_all(a.seed)
    checkpoint = torch.load(a.model_path, map_location="cpu", weights_only=True)
    config = checkpoint["config"]
    model = TabICL(**config)
    model.load_state_dict(checkpoint["state_dict"])
    if not a.unfreeze_encoder:
        for module in [model.col_embedder, model.row_interactor]:
            module.requires_grad_(False)
    model.to(a.device)
    # The held-out validation/calibration/test partitions never supply context Y.
    folds = fold_labels(y, groups, a.folds, a.seed+41)
    pd.DataFrame({"eid": ids, "context_exclusion_fold": folds}).to_csv(out/"inner_folds.csv", index=False)
    context = context_indices(np.arange(len(y)), y, w, a.context_size, np.random.default_rng(a.seed+42))
    pd.DataFrame({"eid": ids[context], "Y": y[context], "IPCW": w[context]}).to_csv(out/"context.csv", index=False)
    optimizer = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                                 lr=a.learning_rate, weight_decay=a.weight_decay)
    scaler = torch.amp.GradScaler("cuda", enabled=a.amp and a.device == "cuda")
    history, best, bad, selected_epoch = [], np.inf, 0, -1
    best_file = out/"selected_weights.ckpt"
    last = time.monotonic()
    # Include the unadapted foundation model as epoch zero. Validation may keep it.
    for epoch in range(a.epochs+1):
        start_time = time.monotonic()
        numerator, denominator, steps = 0., 0., 0
        if epoch:
            rng = np.random.default_rng(a.seed+epoch)
            model.train()  # Frozen encoders use train forward to preserve dtypes.
            for fold in np.unique(folds):
                pool = np.flatnonzero(folds != fold)
                ctx = context_indices(pool, y, w, a.context_size, rng)
                query = rng.permutation(np.flatnonzero(folds == fold))
                gen = EnsembleGenerator(classification=True, n_estimators=1, norm_methods=["none"],
                    feat_shuffle_method="none", class_shuffle_method="none", random_state=a.seed)
                gen.fit(x[ctx], y[ctx])
                for begin in range(0, len(query), a.batch_size):
                    ix = query[begin:begin+a.batch_size]
                    xx, yy = next(iter(gen.transform(x[ix], mode="both").values()))
                    xx = torch.as_tensor(xx, dtype=torch.float32, device=a.device)
                    yy = torch.as_tensor(yy, dtype=torch.float32, device=a.device)
                    target = torch.as_tensor(y[ix], dtype=torch.long, device=a.device)
                    weight = torch.as_tensor(w[ix], dtype=torch.float32, device=a.device)
                    optimizer.zero_grad(set_to_none=True)
                    with torch.autocast(device_type=a.device, dtype=torch.float16,
                                        enabled=a.amp and a.device == "cuda"):
                        logits = model(xx, yy)[0, :, :2]
                        batch_loss = (F.cross_entropy(logits.float(), target, reduction="none")*weight).sum()/weight.sum()
                    if not torch.isfinite(batch_loss):
                        raise FloatingPointError("Nonfinite TabICL training loss")
                    scaler.scale(batch_loss).backward()
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.)
                    scaler.step(optimizer)
                    scaler.update()
                    numerator += batch_loss.item()*float(weight.sum())
                    denominator += float(weight.sum())
                    steps += 1
                    if time.monotonic()-last > 45:
                        log("RUNNING", "TF finetune", f"epoch={epoch}; steps={steps}")
                        last = time.monotonic()
        estimator = _predictor(model, config, x[context], y[context], a)
        valid = predict(estimator, xv, a.predict_batch_size)
        score = float(np.average(loss(yv, valid), weights=wv))
        if not np.isfinite(score):
            raise FloatingPointError("Nonfinite TF validation score")
        improved = score < best-1e-6
        if improved:
            best, bad, selected_epoch = score, 0, epoch
            state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            tmp = best_file.with_suffix(".tmp")
            torch.save(dict(config=config, state_dict=state), tmp)
            tmp.replace(best_file)
        else:
            bad += 1
        history.append(dict(epoch=epoch, steps=steps, train_IPCW_logloss=numerator/denominator if denominator else None,
                            tune_IPCW_logloss=score, selected=improved, seconds=time.monotonic()-start_time))
        pd.DataFrame(history).to_csv(out/"learning_curve.csv", index=False)
        log("EPOCH", "TF finetune", f"{epoch}/{a.epochs}; tune_logloss={score:.5f}; steps={steps}; best_epoch={selected_epoch}")
        del estimator
        model.clear_cache()
        if epoch and bad >= a.patience:
            break
    selected = torch.load(best_file, map_location="cpu", weights_only=True)
    model.load_state_dict(selected["state_dict"])
    changes = sum(not torch.equal(v, checkpoint["state_dict"][k]) for k, v in selected["state_dict"].items())
    dump(out/"training_summary.json", dict(pretrained_parameters=sum(p.numel() for p in model.parameters()),
        trainable_parameters=sum(p.numel() for p in model.parameters() if p.requires_grad),
        gradient_steps=sum(r["steps"] for r in history), best_epoch=selected_epoch,
        selected_tensors_changed=changes, loss="IPCW weighted cross entropy", context_rows=len(context),
        family_context_exclusion=groups is not None, pretrained_baseline_epoch=0,
        context_weighting="weighted sampling without replacement; approximate, not exact IPCW attention"))
    estimator = _predictor(model, config, x[context], y[context], a)
    return estimator, selected, context


def restore_predictor(bundle, a):
    from tabicl._model import TabICL
    checkpoint = bundle["transformer"]
    model = TabICL(**checkpoint["config"])
    model.load_state_dict(checkpoint["state_dict"])
    model.to(a.device)
    return _predictor(model, checkpoint["config"], bundle["context_x"], bundle["context_y"], a)
