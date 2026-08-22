#!/usr/bin/env python3
"""Fine-tune PubMedBERT/BiomedBERT for journal placement and paper-score prediction."""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass

def r2_score_np(y: np.ndarray, p: np.ndarray) -> float:
    mask = np.isfinite(y) & np.isfinite(p)
    if mask.sum() < 3:
        return float("nan")
    y = y[mask]
    p = p[mask]
    ss_res = np.sum((y - p) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

def pearson_np(y: np.ndarray, p: np.ndarray) -> float:
    mask = np.isfinite(y) & np.isfinite(p)
    if mask.sum() < 3:
        return float("nan")
    y = y[mask]
    p = p[mask]
    if np.std(y) == 0 or np.std(p) == 0:
        return float("nan")
    return float(np.corrcoef(y, p)[0, 1])

def metrics(y: np.ndarray, p: np.ndarray, prefix: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    names = ["journal", "paper"]
    for i, name in enumerate(names):
        yy = y[:, i]
        pp = p[:, i]
        mask = np.isfinite(yy) & np.isfinite(pp)
        if mask.sum() == 0:
            continue
        out[f"{prefix}_{name}_mae"] = float(np.mean(np.abs(yy[mask] - pp[mask])))
        out[f"{prefix}_{name}_rmse"] = float(np.sqrt(np.mean((yy[mask] - pp[mask]) ** 2)))
        out[f"{prefix}_{name}_r2"] = r2_score_np(yy, pp)
        out[f"{prefix}_{name}_pearson"] = pearson_np(yy, pp)
    return out

def fit_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_specs = [("all", None)]
    for col in ["domain", "paper_score_source"]:
        if col in df.columns:
            group_specs.extend((f"{col}={k}", v.index) for k, v in df.groupby(col, dropna=False))
    for group, idx in group_specs:
        d = df if idx is None else df.loc[idx]
        row = {"group": group, "n": int(len(d))}
        pairs = [
            ("journal", "journal_score_0_100", "pred_expected_journal_score_0_100"),
            ("paper", "paper_score", "pred_paper_score"),
        ]
        for name, ycol, pcol in pairs:
            if ycol not in d.columns or pcol not in d.columns:
                continue
            y = pd.to_numeric(d[ycol], errors="coerce")
            p = pd.to_numeric(d[pcol], errors="coerce")
            mask = y.notna() & p.notna()
            if mask.sum() < 3:
                continue
            yy = y[mask]
            pp = p[mask]
            row[f"{name}_mae"] = float((yy - pp).abs().mean())
            row[f"{name}_rmse"] = float(np.sqrt(((yy - pp) ** 2).mean()))
            row[f"{name}_pearson"] = float(yy.corr(pp, method="pearson"))
            row[f"{name}_spearman"] = float(yy.corr(pp, method="spearman"))
        rows.append(row)
    return pd.DataFrame(rows)

class PFIDataset:
    def __init__(self, df: pd.DataFrame):
        self.texts = df["blinded_text"].fillna("").astype(str).tolist()
        self.labels = df[["journal_score_0_100", "paper_score"]].astype("float32").values / 100.0
        self.meta = df.drop(columns=["blinded_text"], errors="ignore")

    def __len__(self) -> int:
        return len(self.texts)

    def __getitem__(self, idx: int) -> Dict[str, object]:
        return {"text": self.texts[idx], "labels": self.labels[idx]}

def make_loader(df: pd.DataFrame, tokenizer, batch_size: int, max_length: int, shuffle: bool, num_workers: int = 0):
    import torch
    from torch.utils.data import DataLoader

    dataset = PFIDataset(df)

    def collate(batch: List[Dict[str, object]]):
        texts = [b["text"] for b in batch]
        enc = tokenizer(
            texts,
            truncation=True,
            padding=True,
            max_length=max_length,
            return_tensors="pt",
        )
        enc["labels"] = torch.tensor(np.stack([b["labels"] for b in batch]), dtype=torch.float32)
        return enc

    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, collate_fn=collate, num_workers=num_workers)

def run_eval(model, loader, device):
    import torch

    model.eval()
    ys = []
    ps = []
    with torch.no_grad():
        for batch in loader:
            labels = batch.pop("labels").to(device)
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(**batch)
            pred = torch.clamp(out.logits, 0.0, 1.0)
            ys.append(labels.detach().cpu().numpy() * 100.0)
            ps.append(pred.detach().cpu().numpy() * 100.0)
    if not ys:
        return np.empty((0, 2)), np.empty((0, 2))
    return np.vstack(ys), np.vstack(ps)

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--training", required=True)
    ap.add_argument("--model_name", default="microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext")
    ap.add_argument("--fallback_model_name", default="microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext")
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--metrics_out", required=True)
    ap.add_argument("--valid_pred_out", required=True)
    ap.add_argument("--fit_summary_out", default=None)
    ap.add_argument("--max_length", type=int, default=512)
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--batch_size", type=int, default=16)
    ap.add_argument("--grad_accum", type=int, default=2)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--weight_decay", type=float, default=0.01)
    ap.add_argument("--warmup_ratio", type=float, default=0.06)
    ap.add_argument("--max_train_rows", type=int, default=0)
    ap.add_argument("--num_workers", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fp16", action="store_true")
    args = ap.parse_args()

    set_seed(args.seed)
    Path(args.model_dir).mkdir(parents=True, exist_ok=True)
    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.valid_pred_out).parent.mkdir(parents=True, exist_ok=True)

    try:
        import torch
        from torch.optim import AdamW
        from transformers import AutoConfig, AutoModelForSequenceClassification, AutoTokenizer, get_linear_schedule_with_warmup
    except Exception as e:
        raise RuntimeError(
            "Missing dependencies. Install torch, transformers, accelerate, huggingface_hub, pyarrow, duckdb, pandas."
        ) from e

    df = pd.read_parquet(args.training)
    df = df.dropna(subset=["blinded_text", "journal_score_0_100", "paper_score"])
    if args.max_train_rows and args.max_train_rows > 0 and len(df) > args.max_train_rows:
        df = df.sample(args.max_train_rows, random_state=args.seed)
    if "split" not in df.columns:
        rng = np.random.default_rng(args.seed)
        r = rng.random(len(df))
        df["split"] = np.where(r < 0.1, "test", np.where(r < 0.2, "valid", "train"))

    train_df = df[df["split"].eq("train")].reset_index(drop=True)
    valid_df = df[df["split"].eq("valid")].reset_index(drop=True)
    test_df = df[df["split"].eq("test")].reset_index(drop=True)
    if valid_df.empty:
        valid_df = train_df.sample(min(1000, len(train_df)), random_state=args.seed)
    if test_df.empty:
        test_df = valid_df.copy()

    print(f"Training rows: train={len(train_df):,} valid={len(valid_df):,} test={len(test_df):,}")
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    def load_model(name: str):
        tok = AutoTokenizer.from_pretrained(name)
        cfg = AutoConfig.from_pretrained(name, num_labels=2)
        cfg.problem_type = "regression"
        mod = AutoModelForSequenceClassification.from_pretrained(name, config=cfg)
        return tok, mod

    try:
        tokenizer, model = load_model(args.model_name)
        used_model = args.model_name
    except Exception as first_error:
        print(f"WARNING: failed to load {args.model_name}: {first_error}")
        print(f"Trying fallback model: {args.fallback_model_name}")
        tokenizer, model = load_model(args.fallback_model_name)
        used_model = args.fallback_model_name

    model.to(device)
    train_loader = make_loader(train_df, tokenizer, args.batch_size, args.max_length, True, args.num_workers)
    valid_loader = make_loader(valid_df, tokenizer, args.batch_size, args.max_length, False, args.num_workers)
    test_loader = make_loader(test_df, tokenizer, args.batch_size, args.max_length, False, args.num_workers)

    optim = AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    total_update_steps = math.ceil(len(train_loader) / max(args.grad_accum, 1)) * max(args.epochs, 1)
    warmup_steps = int(total_update_steps * args.warmup_ratio)
    sched = get_linear_schedule_with_warmup(optim, num_warmup_steps=warmup_steps, num_training_steps=total_update_steps)
    scaler = torch.cuda.amp.GradScaler(enabled=bool(args.fp16 and device.type == "cuda"))

    best_valid = float("inf")
    best_state = None
    global_step = 0
    train_history = []
    for epoch in range(1, args.epochs + 1):
        model.train()
        losses = []
        optim.zero_grad(set_to_none=True)
        for step, batch in enumerate(train_loader, 1):
            labels = batch.pop("labels").to(device)
            batch = {k: v.to(device) for k, v in batch.items()}
            with torch.cuda.amp.autocast(enabled=bool(args.fp16 and device.type == "cuda")):
                out = model(**batch, labels=labels)
                loss = out.loss / max(args.grad_accum, 1)
            scaler.scale(loss).backward()
            losses.append(float(loss.detach().cpu()) * max(args.grad_accum, 1))
            if step % args.grad_accum == 0 or step == len(train_loader):
                scaler.unscale_(optim)
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                scaler.step(optim)
                scaler.update()
                sched.step()
                optim.zero_grad(set_to_none=True)
                global_step += 1
        yv, pv = run_eval(model, valid_loader, device)
        valid_mae = float(np.mean(np.abs(yv - pv))) if len(yv) else float("inf")
        train_loss = float(np.mean(losses)) if losses else float("nan")
        print(f"epoch={epoch} train_loss={train_loss:.5f} valid_mean_mae={valid_mae:.3f}")
        train_history.append({"epoch": epoch, "train_loss": train_loss, "valid_mean_mae": valid_mae})
        if valid_mae < best_valid:
            best_valid = valid_mae
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    if best_state is not None:
        model.load_state_dict(best_state)
        model.to(device)

    yv, pv = run_eval(model, valid_loader, device)
    yt, pt = run_eval(model, test_loader, device)
    met = {
        "model_name": used_model,
        "n_train": int(len(train_df)),
        "n_valid": int(len(valid_df)),
        "n_test": int(len(test_df)),
        "max_length": int(args.max_length),
        "epochs": int(args.epochs),
        "batch_size": int(args.batch_size),
        "grad_accum": int(args.grad_accum),
        "lr": float(args.lr),
        "seed": int(args.seed),
        "history": train_history,
    }
    met.update(metrics(yv, pv, "valid"))
    met.update(metrics(yt, pt, "test"))

    model.save_pretrained(args.model_dir)
    tokenizer.save_pretrained(args.model_dir)
    config = {
        "base_model": used_model,
        "targets": ["expected_journal_score_0_100", "predicted_paper_score"],
        "training_file": str(args.training),
        "scale": "model logits are trained on labels divided by 100; prediction script rescales to 0-100",
    }
    Path(args.model_dir, "model_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    Path(args.metrics_out).write_text(json.dumps(met, indent=2), encoding="utf-8")

    valid_out = valid_df.drop(columns=["blinded_text"], errors="ignore").copy()
    valid_out["pred_expected_journal_score_0_100"] = pv[:, 0] if len(pv) else np.nan
    valid_out["pred_paper_score"] = pv[:, 1] if len(pv) else np.nan
    valid_out.to_csv(args.valid_pred_out, index=False)
    if args.fit_summary_out:
        Path(args.fit_summary_out).parent.mkdir(parents=True, exist_ok=True)
        fit_summary_table(valid_out).to_csv(args.fit_summary_out, index=False)
    print(f"Saved model: {args.model_dir}")
    print(f"Saved metrics: {args.metrics_out}")
    print(f"Saved validation predictions: {args.valid_pred_out}")
    if args.fit_summary_out:
        print(f"Saved fit summary: {args.fit_summary_out}")

if __name__ == "__main__":
    main()
