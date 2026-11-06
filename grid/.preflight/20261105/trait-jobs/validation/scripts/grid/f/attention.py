"""GRID's module-token Transformer and selective donor attention.

Adapted from LE8 c1's ContextStack / RowEncoder / QK borrowing design for
continuous traits and baseline binary disease. Values are CSx OOF residuals,
not survival labels; genetic block distances supply the candidate caliper.
Only build fits weights/preprocessing; tune_model selects the checkpoint.
The caller owns the independent gate, calibration audit and outer test split.
"""
from __future__ import annotations

import math
import importlib.util
import os
import sys
from pathlib import Path
import numpy as np

# Required before CUDA creates a cuBLAS handle. Scientific reruns (including
# the test-label isolation audit) must not drift with atomic reduction order.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import torch
from torch import nn
from torch.nn import functional as F

# Shared kernels live in 0f; GRID's encoder and training policy stay in this file.
# Keep one module identity so both projects and saved models use the same classes.
_SHARED_ABM_PATH = Path(__file__).resolve().parents[2] / "0f" / "abm.py"
shared_abm = sys.modules.get("shared_abm")
if shared_abm is not None:
    if Path(getattr(shared_abm, "__file__", "")).resolve() != _SHARED_ABM_PATH:
        raise ImportError(f"shared_abm is already loaded from a different path; expected {_SHARED_ABM_PATH}")
else:
    _spec = importlib.util.spec_from_file_location("shared_abm", _SHARED_ABM_PATH)
    shared_abm = importlib.util.module_from_spec(_spec)
    sys.modules["shared_abm"] = shared_abm
    try:
        _spec.loader.exec_module(shared_abm)
    except BaseException:
        sys.modules.pop("shared_abm", None)
        raise


def execution_device(request):
    torch.use_deterministic_algorithms(True)
    device = torch.device(request)
    if device.type not in {"cpu", "cuda"}:
        raise ValueError("GRID attention device must be cpu or cuda[:index]")
    if device.type == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError(f"CUDA unavailable in {sys.executable}; GRID refuses CPU fallback")
        if device.index is None:
            device = torch.device("cuda", torch.cuda.current_device())
        torch.cuda.set_device(device)
    return device


class ContextBlock(shared_abm.ContextBlock):
    """Preserve GRID's tensor-only API and saved state keys."""
    def forward(self, x):
        return super().forward(x)[0]


class RowEncoder(nn.Module):
    def __init__(self, membership, width=64, heads=4, layers=2, dropout=.1):
        super().__init__()
        self.width, self.heads = width, heads
        self.n_features = len(membership)
        self.n_tokens = max(membership) + 1
        self.register_buffer("membership", torch.tensor(membership, dtype=torch.long))
        self.register_buffer("counts", torch.tensor(np.bincount(membership), dtype=torch.float32))
        self.value_embedding = nn.Parameter(torch.randn(self.n_features, width) * .05)
        self.missing_embedding = nn.Parameter(torch.randn(self.n_features, width) * .02)
        self.token_embedding = nn.Parameter(torch.randn(self.n_tokens, width) * .02)
        self.cls = nn.Parameter(torch.randn(1, 1, width) * .02)
        self.input_norm = nn.LayerNorm(width)
        self.context = nn.Sequential(*[ContextBlock(width, heads, dropout) for _ in range(layers)],
                                     nn.LayerNorm(width))
        self.project = nn.Sequential(nn.Linear(width, width), nn.GELU(), nn.Linear(width, width))
        self.direct = nn.Linear(width, 1)
        self.decoder_weight = nn.Parameter(torch.randn(self.n_features, width) * .05)
        self.decoder_bias = nn.Parameter(torch.zeros(self.n_features))
        self.query = nn.Linear(width, width, bias=False)
        self.key = nn.Linear(width, width, bias=False)
        self.gate = nn.Linear(width, heads)
        self.temperature = nn.Parameter(torch.full((heads,), -1.))
        with torch.no_grad():
            for linear in (self.query, self.key):
                linear.weight.copy_(torch.eye(width) + torch.randn(width, width) * .01)

    def forward(self, x, observed, decode=False):
        x = torch.where(observed, x, torch.zeros_like(x))
        token = x.new_zeros((len(x), self.n_tokens, self.width))
        token.index_add_(1, self.membership, x[:, :, None] * self.value_embedding)
        token.index_add_(1, self.membership, (~observed)[:, :, None] * self.missing_embedding)
        token = self.input_norm(token / self.counts.sqrt()[None, :, None] + self.token_embedding)
        context = self.context(torch.cat([self.cls.expand(len(x), -1, -1), token], dim=1))
        row = context[:, 0]
        z = F.normalize(self.project(row), dim=-1)
        reconstruction = ((context[:, 1:][:, self.membership] * self.decoder_weight).sum(-1)
                          + self.decoder_bias) if decode else None
        return z, self.direct(row).squeeze(-1), reconstruction

    def donor_weights(self, query, donors, distance_ratio, allowed):
        h, d = self.heads, self.width // self.heads
        q = self.query(query * math.sqrt(self.width)).view(-1, h, d)
        k = self.key(donors * math.sqrt(self.width)).view(len(query), -1, h, d)
        temperature = .02 + 1.98 * torch.sigmoid(self.temperature)
        logits = (q[:, None] * k).sum(-1) / (math.sqrt(d) * temperature)
        logits = logits - .5 * distance_ratio[:, :, None].square()
        # Empty neighborhoods produce exactly zero weights and finite gradients.
        head = shared_abm.safe_softmax(logits, allowed[:, :, None], dim=1)
        return (head * torch.softmax(self.gate(query), -1)[:, None]).sum(-1)


def preflight(device="cuda"):
    dev = execution_device(device)
    # Exercise the real encoder, donor Q/K, gate, decoder and backward kernels.
    with torch.random.fork_rng(devices=[torch.cuda.current_device()] if dev.type == "cuda" else []):
        torch.manual_seed(991)
        model = RowEncoder([0, 0, 1, 1, 2, 2], 16, 2, 1, 0).to(dev)
        x = torch.randn(4, 6, device=dev)
        z, direct, decoded = model(x, torch.ones_like(x, dtype=torch.bool), True)
        weights = model.donor_weights(z, torch.stack([z.roll(j, 0) for j in (1, 2, 3)], dim=1),
            torch.zeros(4, 3, device=dev), torch.ones(4, 3, dtype=torch.bool, device=dev))
        loss = direct.square().mean() + decoded.square().mean() + (weights * torch.arange(3, device=dev)).square().mean()
        loss.backward()
        if not torch.isfinite(loss) or not all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters()):
            raise RuntimeError("GRID attention forward/backward preflight failed")
        for name in ("query.weight", "key.weight", "context.0.qkv.weight"):
            grad = dict(model.named_parameters())[name].grad
            if grad is None or not torch.any(grad != 0):
                raise RuntimeError(f"GRID attention preflight has no {name} gradient")
        if dev.type == "cuda":
            torch.cuda.synchronize(dev)
    return dict(python=sys.executable, device=str(dev), torch=torch.__version__,
        cuda_version=torch.version.cuda, gpu=torch.cuda.get_device_name(dev) if dev.type == "cuda" else None,
        architecture="module_token_ContextStack_QK", forward_backward="passed")


def neighbors(*args, **kwargs):
    return shared_abm.neighbors(*args, **kwargs)


class SelectiveAttention:
    """CPU-serializable frozen state; CUDA models are rebuilt lazily on inference."""
    def __getstate__(self):
        return {k: v for k, v in vars(self).items() if k not in {"_model", "_memory"}}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._model = self._memory = None

    def arrays(self, frame):
        raw = frame[self.columns].to_numpy(float)
        observed = np.isfinite(raw)
        x = (np.where(observed, raw, self.median) - self.center) / self.scale
        return np.clip(x, -20, 20).astype(np.float32), observed

    def runtime(self):
        device = execution_device(self.config["device"])
        if getattr(self, "_model", None) is None:
            self._model = RowEncoder(**self.architecture).to(device)
            self._model.load_state_dict({k: torch.as_tensor(v, device=device) for k, v in self.state.items()})
            self._model.eval()
            self._memory = torch.as_tensor(self.memory, device=device)
        return self._model, self._memory, device

    @torch.no_grad()
    def encode(self, model, x, observed, device):
        model.eval()
        output = []
        for start in range(0, len(x), self.config["attention_batch"]):
            xx = torch.as_tensor(x[start:start+self.config["attention_batch"]], device=device)
            mask = torch.as_tensor(observed[start:start+self.config["attention_batch"]], device=device)
            output.append(model(xx, mask)[0])
        return torch.cat(output)

    def fit(self, bank, build, tune, tune_baseline, y_tune, blocks, config):
        self.config = dict(config)
        device = execution_device(config["device"])
        torch.manual_seed(config["seed"])
        if device.type == "cuda":
            torch.cuda.manual_seed_all(config["seed"])
            torch.cuda.reset_peak_memory_stats(device)
        self.columns, membership = [], []
        for cols in blocks.values():
            token = len(set(membership))
            for column in cols:
                if column not in self.columns:
                    self.columns.append(column)
                    membership.append(token)
        raw = build[self.columns].to_numpy(float)
        self.median = np.array([np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else 0. for v in raw.T])
        filled = np.where(np.isfinite(raw), raw, self.median)
        self.center, self.scale = filled.mean(0), np.maximum(filled.std(0), 1e-8)
        x, observed = self.arrays(build)
        vx, vm = self.arrays(tune)
        self.architecture = dict(membership=membership, width=config["attention_width"],
            heads=config["attention_heads"], layers=config["attention_layers"], dropout=config["attention_dropout"])
        model = RowEncoder(**self.architecture).to(device)
        initial = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        optimizer = torch.optim.AdamW(model.parameters(), lr=config["attention_lr"], weight_decay=1e-4)
        k = min(max(config["k_grid"]), len(build))
        j, d, _ = neighbors(bank.z, bank.z, bank.ids, bank.ids, bank.groups, bank.groups,
            bank.tie, k, device, config["attention_batch"], config["donor_block"],
            bank_folds=bank.donor_folds, query_folds=bank.donor_folds)
        if np.any((j >= 0) & (bank.donor_folds[np.maximum(j, 0)] == bank.donor_folds[:, None])):
            raise AssertionError("An attention training donor crossed its excluded fold")
        vj, vd = bank.query(bank.geometry.transform(tune), tune.eid.to_numpy(str), tune.family_id.to_numpy(str), k)
        radius = max(config["radius_grid"])
        # Standardized residual objectives give height and LDL comparable scales.
        scale = max(float(np.std(bank.residual)), 1e-6)
        target = bank.residual / scale
        vt = (np.asarray(y_tune) - np.asarray(tune_baseline)) / scale
        self.history, steps = [], 0
        named_parameters = list(model.named_parameters())
        gradient_audit = shared_abm.GradientTracker(named_parameters, device)
        rng = np.random.default_rng(config["seed"] + 1901)
        best, patience = None, 0

        _, family_codes = np.unique(bank.groups, return_inverse=True)
        def borrow(query, memory, indices, dist, query_observed):
            allowed = (indices >= 0) & np.isfinite(dist) & (dist <= bank.radius * radius) & (bank.radius > 0)
            allowed &= (1-query_observed.mean(1) <= config["max_missing_fraction"])[:, None]
            safe = np.maximum(indices, 0)
            ratio = np.where(allowed, dist / max(bank.radius, 1e-12), 0.)
            weights = model.donor_weights(query, memory[torch.as_tensor(safe, device=device)],
                torch.as_tensor(ratio, dtype=torch.float32, device=device), torch.as_tensor(allowed, device=device))
            labels = torch.as_tensor(target[safe], dtype=torch.float32, device=device)
            # Family ESS and geometric support use exactly the inference formula.
            codes = torch.as_tensor(family_codes[safe], device=device)
            ess = shared_abm.family_effective_size(weights, codes)
            nearest = np.min(np.where(allowed, dist, np.inf), axis=1)
            support = torch.as_tensor(np.exp(-.5*(nearest/max(bank.radius, 1e-12))**2), dtype=torch.float32, device=device)
            enough = torch.as_tensor(allowed.sum(1) >= config["min_matches"], device=device) & (ess >= config["min_ess"])
            fraction = support * ess / (ess + config["prior_strength"]).clamp_min(1e-20)
            return (weights * labels).sum(1) * fraction * enough

        for epoch in range(1, config["attention_epochs"] + 1):
            memory = self.encode(model, x, observed, device).detach()
            model.train()
            losses = []
            for ix in np.array_split(rng.permutation(len(x)), max(1, math.ceil(len(x)/config["attention_batch"]))):
                xx = torch.as_tensor(x[ix], device=device)
                mask = torch.as_tensor(observed[ix], device=device)
                z, direct, _ = model(xx, mask)
                pred = borrow(z, memory, j[ix], d[ix], observed[ix])
                yy = torch.as_tensor(target[ix], dtype=torch.float32, device=device)
                loss = F.mse_loss(pred, yy) + .25 * F.mse_loss(direct, yy)
                if config["attention_reconstruction"]:
                    hidden = (torch.rand_like(xx) < .15) & mask
                    _, _, decoded = model(xx, mask & ~hidden, True)
                    reconstruction = shared_abm.masked_mse(decoded, xx, hidden)
                    loss = loss + config["attention_reconstruction"] * reconstruction
                if not torch.isfinite(loss):
                    raise RuntimeError("Nonfinite GRID attention training loss")
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1., error_if_nonfinite=True)
                gradient_audit.update()
                optimizer.step()
                losses.append(loss.detach())
                steps += 1
            gradient_audit.flush()
            memory = self.encode(model, x, observed, device).detach()
            vz = self.encode(model, vx, vm, device)
            with torch.no_grad():
                predicted = []
                for start in range(0, len(vx), config["attention_batch"]):
                    stop = start + config["attention_batch"]
                    predicted.append(borrow(vz[start:stop], memory, vj[start:stop], vd[start:stop], vm[start:stop]).cpu().numpy())
                loss = float(np.mean((np.concatenate(predicted) - vt)**2))
            row = dict(epoch=epoch, train_loss=float(torch.stack(losses).mean()), validation_loss=loss)
            self.history.append(row)
            if config.get("verbose", True):
                print(f"[GRID] attention epoch={epoch} device={device} validation_loss={loss:.6g} steps={steps}", flush=True)
            if best is None or loss < best[0] - 1e-9:
                best = (loss, {k: v.detach().cpu().numpy().copy() for k, v in model.state_dict().items()}, epoch)
                patience = 0
            else:
                patience += 1
                if patience >= config["attention_patience"]:
                    break
        self.state, selected_epoch = best[1:]
        model.load_state_dict({k: torch.as_tensor(v, device=device) for k, v in self.state.items()})
        self.memory = self.encode(model, x, observed, device).cpu().numpy()
        self._model = self._memory = None
        self.fit_status = dict(architecture="module_token_ContextStack_QK", actual_device=str(device),
            precision="fp32", gradient_steps=steps, selected_epoch=selected_epoch,
            selected_tensor_changed_count=sum(not np.array_equal(initial[k].numpy(), v) for k, v in self.state.items()),
            trainable_parameters=sum(p.numel() for p in model.parameters()),
            gradient_names=sorted(gradient_audit.flush()),
            train_n=len(build), tune_n=len(tune), donor_fold_exclusion=True,
            gpu_peak_allocated_bytes=int(torch.cuda.max_memory_allocated(device)) if device.type == "cuda" else 0,
            value_source="family_out_of_fold_CSx_residual", checkpoint_metric="tune_model_standardized_residual_MSE")
        return self

    @torch.no_grad()
    def weights(self, frame, indices, distances, valid, radius):
        model, memory, device = self.runtime()
        x, mask = self.arrays(frame)
        weights = np.zeros(indices.shape, np.float64)
        batch = self.config["attention_batch"]
        for start in range(0, len(frame), batch):
            stop = start + batch
            z = model(torch.as_tensor(x[start:stop], device=device),
                      torch.as_tensor(mask[start:stop], device=device))[0]
            allowed = valid[start:stop]
            ratios = np.where(allowed, distances[start:stop] / max(radius, 1e-12), 0.)
            donors = memory[torch.as_tensor(np.maximum(indices[start:stop], 0), device=device)]
            weights[start:stop] = model.donor_weights(z, donors,
                torch.as_tensor(ratios, dtype=torch.float32, device=device),
                torch.as_tensor(allowed, device=device)).cpu().numpy()
        # Float64 normalization preserves exact donor decomposition downstream.
        total = weights.sum(1, keepdims=True)
        return np.divide(weights, total, out=np.zeros_like(weights), where=total > 0)
