"""Continuous-time latent health model, implemented independently in PyTorch.

Reference: HealthFlux v3, doi:10.64898/2026.09.19.26363460.
This is an independently specified implementation of the published state-space
idea, not the authors' code or weights. Numerical integration uses fixed-step
RK4. Numeric MRI-derived features are supported as an ordinary modality; this
module does not silently substitute them for the paper's pretrained raw-3D-MRI
encoder. The explicit implementation choices are recorded in model_cfg.

All times are attained ages in YEARS. next_dt is an elapsed interval in YEARS.
Token zero means no clinical token/padding and is never a clinical target.
An observation mask is not inferred from a zero numeric measurement.
"""
from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import torch
from torch import Tensor, nn
from torch.nn import functional as F


DEFAULTS = {
    "deterministic_dim": 256,
    "stochastic_dim": 128,
    "embedding_dim": 256,
    "context_hidden_dim": 128,
    "age_dim": 64,
    "hidden_dim": 384,
    "dropout": 0.1,
    "max_ode_step_years": 1.0,
    "max_ode_steps": 256,
    "kl_weight": 0.01,
    "kl_balance": 0.8,
    "free_nats": 0.0,
    "time_weight": 1.0,
    "no_ode": False,
    "no_update": False,
    "age_scale": 100.0,
}


def _mlp(inp: int, hidden: int, out: int, dropout: float = 0.0) -> nn.Sequential:
    return nn.Sequential(nn.Linear(inp, hidden), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden, out))


class AgeEncoder(nn.Module):
    """Smooth attained-age basis; no dataset-derived normalization is fitted."""

    def __init__(self, out_dim: int, scale: float = 100.0):
        super().__init__()
        self.scale = float(scale)
        self.register_buffer("frequencies", 2.0 ** torch.arange(6, dtype=torch.float32))
        self.net = nn.Sequential(nn.Linear(15, out_dim), nn.SiLU(), nn.Linear(out_dim, out_dim))

    def forward(self, age: Tensor) -> Tensor:
        a = age / self.scale
        phase = a.unsqueeze(-1) * self.frequencies * (2 * math.pi)
        basis = torch.cat((a.unsqueeze(-1), a.square().unsqueeze(-1),
                           torch.log1p(age.clamp_min(0)).unsqueeze(-1) / math.log1p(self.scale),
                           phase.sin(), phase.cos()), dim=-1)
        return self.net(basis)


class FluxModel(nn.Module):
    """ODE deterministic state + Gaussian stochastic belief + observation GRU.

    Constructor
    -----------
    model_cfg: keys listed in DEFAULTS; unknown keys are rejected.
    n_tokens: entire vocabulary size INCLUDING token 0 and unknown token 1.
    context_dim: raw permanent-context dimension (before value/mask encoding).
    modality_dims: ordered mapping from modality names to numeric feature counts.

    Batch
    -----
    times[B,T], valid[B,T], tokens[B,T], context[B,C], context_mask[B,C]
    numeric[name] = {values[B,T,F_name], mask[B,T,F_name]}
    Optional training labels: next_event[B,T], next_dt[B,T], label_mask[B,T].
    Observations must be ordered in non-decreasing time with right padding.
    Same-day diagnosis ordering is not an observed temporal ordering: labels
    with non-positive next_dt are excluded from the next-event objective.
    """

    def __init__(self, model_cfg: Mapping[str, Any], n_tokens: int,
                 context_dim: int, modality_dims: Mapping[str, int]):
        super().__init__()
        unknown = set(model_cfg) - set(DEFAULTS)
        if unknown:
            raise ValueError("Unknown model configuration keys: " + ", ".join(sorted(unknown)))
        self.model_cfg = {**DEFAULTS, **dict(model_cfg)}
        cfg = self.model_cfg
        for key in ("deterministic_dim", "stochastic_dim", "embedding_dim", "context_hidden_dim", "age_dim", "hidden_dim", "max_ode_steps"):
            if int(cfg[key]) != cfg[key] or int(cfg[key]) < 1:
                raise ValueError(f"model.{key} must be a positive integer")
        for key in ("no_ode", "no_update"):
            if not isinstance(cfg[key], bool):
                raise ValueError(f"model.{key} must be a YAML boolean")
        if not 0 <= float(cfg["dropout"]) < 1:
            raise ValueError("model.dropout must be in [0,1)")
        if not 0 <= float(cfg["kl_balance"]) <= 1:
            raise ValueError("model.kl_balance must be in [0,1]")
        for key in ("max_ode_step_years", "age_scale"):
            if not math.isfinite(float(cfg[key])) or float(cfg[key]) <= 0:
                raise ValueError(f"model.{key} must be finite and positive")
        for key in ("kl_weight", "free_nats", "time_weight"):
            if not math.isfinite(float(cfg[key])) or float(cfg[key]) < 0:
                raise ValueError(f"model.{key} must be finite and non-negative")
        self.n_tokens = int(n_tokens)
        self.context_dim = int(context_dim)
        if self.n_tokens < 2 or self.context_dim < 0:
            raise ValueError("n_tokens must include PAD and UNK, and context_dim cannot be negative")
        self.modality_dims = {str(k): int(v) for k, v in modality_dims.items()}
        if any(v < 1 for v in self.modality_dims.values()):
            raise ValueError("All enabled numeric modalities need at least one feature")
        self.modality_names = tuple(self.modality_dims)
        self.h_dim = int(cfg["deterministic_dim"])
        self.z_dim = int(cfg["stochastic_dim"])
        self.state_dim = self.h_dim + self.z_dim
        self.embedding_dim = int(cfg["embedding_dim"])
        self.c_dim = int(cfg["context_hidden_dim"])
        self.a_dim = int(cfg["age_dim"])
        hidden = int(cfg["hidden_dim"])
        drop = float(cfg["dropout"])
        self.context_encoder = _mlp(2 * self.context_dim, self.c_dim, self.c_dim) if self.context_dim else None
        self.empty_context = nn.Parameter(torch.zeros(self.c_dim)) if not self.context_dim else None
        self.age_encoder = AgeEncoder(self.a_dim, float(cfg["age_scale"]))
        self.initial_h = _mlp(self.c_dim + self.a_dim, hidden, self.h_dim)
        self.token_embedding = nn.Embedding(self.n_tokens, self.embedding_dim, padding_idx=0)
        # ModuleList permits scientific modality names containing dots.
        self.numeric_encoders = nn.ModuleList([
            _mlp(2 * self.modality_dims[name], self.embedding_dim, self.embedding_dim, drop)
            for name in self.modality_names
        ])
        self.modality_gates = nn.Parameter(torch.zeros(len(self.modality_names) + 1))
        self.observation_norm = nn.LayerNorm(self.embedding_dim)
        dynamic_input = self.h_dim + self.c_dim + self.a_dim
        # No dropout in the vector field: repeated RK stages evaluate the same ODE.
        self.vector_field = nn.Sequential(nn.LayerNorm(self.h_dim), nn.Linear(self.h_dim, hidden))
        self.vector_context = nn.Linear(self.c_dim + self.a_dim, hidden, bias=False)
        self.vector_output = nn.Sequential(nn.SiLU(), nn.Linear(hidden, self.h_dim), nn.Tanh())
        nn.init.normal_(self.vector_output[1].weight, std=0.01)
        nn.init.zeros_(self.vector_output[1].bias)
        self.prior_net = _mlp(dynamic_input, hidden, 2 * self.z_dim)
        posterior_input = dynamic_input + self.embedding_dim
        self.posterior_net = _mlp(posterior_input, hidden, 2 * self.z_dim, drop)
        self.posterior_gate = nn.Linear(posterior_input, self.z_dim)
        self.jump = nn.GRUCell(self.embedding_dim + self.z_dim + self.c_dim + self.a_dim, self.h_dim)
        self.state_norm = nn.LayerNorm(self.state_dim)
        self.event_head = nn.Linear(self.state_dim, self.n_tokens)
        self.time_head = _mlp(self.state_dim, hidden, 1, drop)

    @property
    def device(self) -> torch.device:
        return self.token_embedding.weight.device

    @property
    def dtype(self) -> torch.dtype:
        return self.token_embedding.weight.dtype

    def encode_context(self, values: Tensor, mask: Tensor | None = None) -> Tensor:
        if values.ndim != 2 or values.shape[1] != self.context_dim:
            raise ValueError(f"context must have shape [B,{self.context_dim}]")
        values = values.to(device=self.device, dtype=self.dtype)
        if self.context_encoder is None:
            return self.empty_context.unsqueeze(0).expand(values.shape[0], -1)
        observed = torch.isfinite(values) if mask is None else mask.to(device=self.device, dtype=torch.bool) & torch.isfinite(values)
        if observed.shape != values.shape:
            raise ValueError("context_mask must match context")
        clean = torch.where(observed, values, torch.zeros_like(values))
        return self.context_encoder(torch.cat((clean, observed.to(self.dtype)), -1))

    @staticmethod
    def _normal(params: Tensor) -> tuple[Tensor, Tensor]:
        mean, log_std = params.chunk(2, dim=-1)
        return mean, log_std.clamp(-5.0, 2.0)

    def _prior(self, h: Tensor, context: Tensor, age: Tensor) -> tuple[Tensor, Tensor]:
        return self._normal(self.prior_net(torch.cat((h, context, self.age_encoder(age)), -1)))

    @staticmethod
    def _sample(mean: Tensor, log_std: Tensor, sample: bool) -> Tensor:
        return mean + log_std.exp() * torch.randn_like(mean) if sample else mean

    def initial_state(self, context: Tensor, age: Tensor, context_mask: Tensor | None = None,
                      sample: bool = False) -> dict[str, Tensor]:
        encoded = self.encode_context(context, context_mask)
        age = torch.as_tensor(age, device=self.device, dtype=self.dtype)
        if age.ndim == 0:
            age = age.expand(context.shape[0])
        if age.shape != (context.shape[0],) or not bool(torch.isfinite(age).all()) or bool((age < 0).any()):
            raise ValueError("Initial age must be a finite non-negative scalar or [B] vector")
        h = self.initial_h(torch.cat((encoded, self.age_encoder(age)), -1))
        mu, log_std = self._prior(h, encoded, age)
        return {"h": h, "z": self._sample(mu, log_std, sample), "age": age,
                "context": encoded, "initialized": torch.ones_like(age, dtype=torch.bool)}

    def _drift(self, h: Tensor, age: Tensor, context: Tensor) -> Tensor:
        return self.vector_output(self.vector_field(h) + self.vector_context(torch.cat((context, self.age_encoder(age)), -1)))

    def _integrate(self, h: Tensor, age: Tensor, elapsed: Tensor, context: Tensor) -> Tensor:
        if self.model_cfg["no_ode"]:
            return h
        max_dt = float(elapsed.detach().max().cpu()) if elapsed.numel() else 0.0
        steps = int(math.ceil(max_dt / float(self.model_cfg["max_ode_step_years"])))
        if steps == 0:
            return h
        if steps > int(self.model_cfg["max_ode_steps"]):
            raise ValueError(f"ODE interval needs {steps} RK4 steps; verify dates or raise model.max_ode_steps")
        dt = elapsed / steps
        step = dt.unsqueeze(-1)
        current = age
        result = h
        for _ in range(steps):
            k1 = self._drift(result, current, context)
            k2 = self._drift(result + step * k1 / 2, current + dt / 2, context)
            k3 = self._drift(result + step * k2 / 2, current + dt / 2, context)
            k4 = self._drift(result + step * k3, current + dt, context)
            result = result + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            current = current + dt
        return result

    def advance_state(self, state: Mapping[str, Tensor], elapsed: Tensor | float,
                      sample: bool = False, *, age: Tensor | None = None,
                      context: Tensor | None = None) -> dict[str, Tensor]:
        """No-observation simulation. context, if supplied, is already encoded.

        The state is returned as a new dictionary. A zero interval is exactly
        identity, including the stochastic component. Negative intervals fail
        instead of silently incorporating a future record into a past query.
        """
        current_age = state["age"] if age is None else age
        encoded = state["context"] if context is None else context
        elapsed = torch.as_tensor(elapsed, device=state["h"].device, dtype=state["h"].dtype)
        if elapsed.ndim == 0:
            elapsed = elapsed.expand_as(current_age)
        if elapsed.shape != current_age.shape or not bool(torch.isfinite(elapsed).all()) or bool((elapsed < -1e-5).any()):
            raise ValueError("elapsed must be finite, non-negative and scalar or [B]")
        elapsed = elapsed.clamp_min(0)
        h = self._integrate(state["h"], current_age, elapsed, encoded)
        new_age = current_age + elapsed
        mu, log_std = self._prior(h, encoded, new_age)
        proposal = self._sample(mu, log_std, sample)
        z = torch.where((elapsed > 0).unsqueeze(-1), proposal, state["z"])
        return {"h": h, "z": z, "age": new_age, "context": encoded,
                "initialized": state.get("initialized", torch.ones_like(new_age, dtype=torch.bool))}

    @staticmethod
    def _kl(q_mu: Tensor, q_log_std: Tensor, p_mu: Tensor, p_log_std: Tensor) -> Tensor:
        # Compute in float32 even under mixed precision.
        qm, ql, pm, pl = (x.float() for x in (q_mu, q_log_std, p_mu, p_log_std))
        return (pl - ql + 0.5 * ((2 * (ql - pl)).exp() + (qm - pm).square() * (-2 * pl).exp() - 1)).sum(-1).clamp_min(0)

    def _observe(self, state: Mapping[str, Tensor], observation: Tensor,
                 present: Tensor, sample: bool) -> tuple[dict[str, Tensor], Tensor, Tensor]:
        h, context, age = state["h"], state["context"], state["age"]
        age_emb = self.age_encoder(age)
        prior_mu, prior_log_std = self._prior(h, context, age)
        posterior_input = torch.cat((h, context, age_emb, observation), -1)
        candidate_mu, candidate_log_std = self._normal(self.posterior_net(posterior_input))
        gate = torch.sigmoid(self.posterior_gate(posterior_input))
        posterior_mu = (1 - gate) * prior_mu + gate * candidate_mu
        posterior_log_std = torch.log((1 - gate) * prior_log_std.exp() + gate * candidate_log_std.exp())
        proposed_z = self._sample(posterior_mu, posterior_log_std, sample)
        # The v3 no-update ablation removes recurrent revision, not observation
        # conditioning: the instantaneous stochastic posterior and KL remain.
        proposed_h = h if self.model_cfg["no_update"] else self.jump(
            torch.cat((observation, proposed_z, context, age_emb), -1), h)
        has = present.unsqueeze(-1)
        new = {**state, "h": torch.where(has, proposed_h, h),
               "z": torch.where(has, proposed_z, state["z"])}
        raw = self._kl(posterior_mu, posterior_log_std, prior_mu, prior_log_std)
        a = float(self.model_cfg["kl_balance"])
        balanced = (a * self._kl(posterior_mu.detach(), posterior_log_std.detach(), prior_mu, prior_log_std)
                    + (1 - a) * self._kl(posterior_mu, posterior_log_std, prior_mu.detach(), prior_log_std.detach()))
        balanced = balanced.clamp_min(float(self.model_cfg["free_nats"]))
        return new, raw * present, balanced * present

    def state_features(self, state: Mapping[str, Tensor]) -> Tensor:
        return self.state_norm(torch.cat((state["h"], state["z"]), -1))

    def predict_next(self, state: Mapping[str, Tensor]) -> tuple[Tensor, Tensor]:
        """Return clinical logits and predicted log(1 + interval in years)."""
        features = self.state_features(state)
        return self.event_head(features), F.softplus(self.time_head(features).squeeze(-1))

    def observe_event(self, state: Mapping[str, Tensor], tokens: Tensor,
                      sample: bool = False) -> dict[str, Tensor]:
        """Insert generated clinical events AFTER advance_state in a rollout."""
        tokens = torch.as_tensor(tokens, device=self.device, dtype=torch.long)
        if tokens.shape != state["age"].shape or bool(((tokens < 0) | (tokens >= self.n_tokens)).any()):
            raise ValueError("Event tokens must be an in-vocabulary [B] vector")
        observation = self.observation_norm(self.token_embedding(tokens) * torch.sigmoid(self.modality_gates[0]))
        result, _, _ = self._observe(state, observation, tokens > 0, sample)
        return result

    def _observations(self, batch: Mapping[str, Any], tokens: Tensor, valid: Tensor) -> tuple[Tensor, Tensor]:
        numeric = batch.get("numeric", {})
        unknown = set(numeric) - set(self.modality_names)
        if unknown:
            raise ValueError("Unconfigured numeric modalities: " + ", ".join(sorted(unknown)))
        present = (tokens > 0) & valid
        count = present.to(self.dtype)
        fusion = self.token_embedding(tokens) * torch.sigmoid(self.modality_gates[0])
        for index, name in enumerate(self.modality_names):
            if name not in numeric:
                continue
            values = numeric[name]["values"].to(device=self.device, dtype=self.dtype)
            mask = numeric[name]["mask"].to(device=self.device, dtype=torch.bool)
            expected = (*tokens.shape, self.modality_dims[name])
            if values.shape != expected or mask.shape != expected:
                raise ValueError(f"numeric[{name}] values/mask must have shape {expected}")
            mask = mask & torch.isfinite(values) & valid.unsqueeze(-1)
            available = mask.any(-1)
            # Avoid applying a large omics encoder to EHR-only observations.
            if bool(available.any()):
                clean = torch.where(mask, values, torch.zeros_like(values))
                encoded = self.numeric_encoders[index](torch.cat((clean[available], mask[available].to(self.dtype)), -1))
                addition = torch.zeros_like(fusion)
                addition[available] = encoded * torch.sigmoid(self.modality_gates[index + 1])
                fusion = fusion + addition
                count = count + available.to(self.dtype)
                present = present | available
        fusion = fusion / count.clamp_min(1).sqrt().unsqueeze(-1)
        return self.observation_norm(fusion), present

    def forward(self, batch: Mapping[str, Any], sample: bool | None = None,
                initial_state: Mapping[str, Tensor] | None = None,
                return_sequence: bool = True, compute_loss: bool = True) -> dict[str, Any]:
        sample = self.training if sample is None else bool(sample)
        times = batch["times"].to(device=self.device, dtype=self.dtype)
        valid = batch["valid"].to(device=self.device, dtype=torch.bool)
        tokens = batch["tokens"].to(device=self.device, dtype=torch.long)
        if times.ndim != 2 or times.shape != valid.shape or times.shape != tokens.shape:
            raise ValueError("times, valid and tokens must have matching [B,T] shapes")
        if bool((valid[:, 1:] & ~valid[:, :-1]).any()):
            raise ValueError("Sequences require contiguous valid nodes followed by right padding")
        if bool((~torch.isfinite(times[valid])).any()) or bool((times[valid] < 0).any()):
            raise ValueError("Observed ages must be finite and non-negative")
        if bool(((tokens < 0) | (tokens >= self.n_tokens)).any()):
            raise ValueError("Clinical tokens outside vocabulary")
        if times.shape[1] > 1 and bool(((times[:, 1:] < times[:, :-1] - 1e-5) & valid[:, 1:]).any()):
            raise ValueError("Observation times decrease within a person")
        if initial_state is None:
            if "start_age" in batch:
                first_age = torch.as_tensor(batch["start_age"], device=self.device, dtype=self.dtype)
            elif times.shape[1] and bool(valid.any(-1).all()):
                first_age = times[:, 0]
            else:
                raise ValueError("Empty histories require explicit start_age for every person")
            state = self.initial_state(batch["context"], first_age, batch.get("context_mask"), sample=sample)
        else:
            state = dict(initial_state)
            if state["h"].shape != (times.shape[0], self.h_dim):
                raise ValueError("Initial recurrent state has incompatible batch/dimension")
        observation, present = self._observations(batch, tokens, valid)
        features = []
        raw_kl = []
        balanced_kl = []
        keep_sequence = return_sequence or compute_loss
        for k in range(times.shape[1]):
            elapsed = torch.where(valid[:, k], times[:, k] - state["age"], torch.zeros_like(state["age"]))
            state = self.advance_state(state, elapsed, sample=sample)
            state, kl, bkl = self._observe(state, observation[:, k], present[:, k] & valid[:, k], sample)
            raw_kl.append(kl)
            balanced_kl.append(bkl)
            if keep_sequence:
                features.append(self.state_features(state))
        final = self.state_features(state)
        feature_tensor = torch.stack(features, 1) if features else final.new_empty((times.shape[0], 0, self.state_dim))
        logits = self.event_head(feature_tensor) if keep_sequence else None
        time_prediction = F.softplus(self.time_head(feature_tensor).squeeze(-1)) if keep_sequence else None
        zero = final.sum() * 0.0
        if keep_sequence:
            # A rank can have no eligible next-event targets in this minibatch.
            # Keep its predictive heads connected with zero gradients for DDP.
            zero = zero + logits.sum() * 0.0 + time_prediction.sum() * 0.0
        raw_sum = torch.stack(raw_kl, 1).sum() if raw_kl else zero
        balanced_sum = torch.stack(balanced_kl, 1).sum() if balanced_kl else zero
        kl_count = present.sum()
        event_sum = zero
        time_sum = zero
        event_count = valid.new_tensor(0, dtype=torch.long)
        effective_mask = torch.zeros_like(valid)
        if compute_loss and "label_mask" in batch:
            target = batch["next_event"].to(device=self.device, dtype=torch.long)
            delta = batch["next_dt"].to(device=self.device, dtype=self.dtype)
            supplied = batch["label_mask"].to(device=self.device, dtype=torch.bool)
            if target.shape != valid.shape or delta.shape != valid.shape or supplied.shape != valid.shape:
                raise ValueError("Next-event labels must have shape [B,T]")
            if bool((supplied & ((target <= 0) | (target >= self.n_tokens) | ~torch.isfinite(delta) | (delta < 0))).any()):
                raise ValueError("Enabled training targets must be real vocabulary tokens with finite non-negative intervals")
            effective_mask = valid & supplied & (delta > 0)
            event_count = effective_mask.sum()
            if bool(event_count):
                # PAD is not a competing clinical outcome, even during training.
                selected = logits[effective_mask].clone()
                selected[:, 0] = torch.finfo(selected.dtype).min
                event_sum = F.cross_entropy(selected, target[effective_mask], reduction="sum")
                time_sum = F.smooth_l1_loss(time_prediction[effective_mask], torch.log1p(delta[effective_mask]), reduction="sum")
        event_mean = event_sum / event_count.clamp_min(1)
        time_mean = time_sum / event_count.clamp_min(1)
        kl_mean = balanced_sum / kl_count.clamp_min(1)
        loss = event_mean + float(self.model_cfg["time_weight"]) * time_mean + float(self.model_cfg["kl_weight"]) * kl_mean
        return {"loss": loss, "loss_event": event_mean, "loss_time": time_mean,
                "loss_kl": raw_sum / kl_count.clamp_min(1), "loss_kl_balanced": kl_mean,
                "event_sum": event_sum, "time_sum": time_sum, "kl_sum": raw_sum,
                "kl_balanced_sum": balanced_sum, "event_count": event_count,
                "kl_count": kl_count, "label_mask": effective_mask,
                "logits": logits, "time_prediction": time_prediction,
                "features": feature_tensor if return_sequence else None,
                "final_features": final, "final_state": state}

    def encode_history(self, batch: Mapping[str, Any], sample: bool = False,
                       chunk_size: int = 96) -> dict[str, Tensor]:
        """Encode every supplied observation, maintaining state across chunks.

        Caller should use model.eval() and torch.no_grad() for frozen inference.
        query_age, if supplied, advances to that attained age AFTER processing
        the prefix. It never truncates a prefix that already contains later data.
        """
        if chunk_size < 1:
            raise ValueError("chunk_size must be positive")
        length = batch["times"].shape[1]
        state = None
        for start in range(0, max(length, 1), chunk_size):
            stop = min(length, start + chunk_size)
            piece = {k: batch[k] for k in ("context", "context_mask", "start_age") if k in batch}
            for key in ("times", "valid", "tokens"):
                piece[key] = batch[key][:, start:stop]
            piece["numeric"] = {name: {key: arr[:, start:stop] for key, arr in value.items()}
                                for name, value in batch.get("numeric", {}).items()}
            # Later chunks can be all padding for a subset; carry state handles it.
            state = self.forward(piece, sample=sample, initial_state=state,
                                 return_sequence=False, compute_loss=False)["final_state"]
        if "query_age" in batch:
            query_age = torch.as_tensor(batch["query_age"], device=self.device, dtype=self.dtype)
            state = self.advance_state(state, query_age - state["age"], sample=sample)
        return state


class RiskHead(nn.Module):
    """Endpoint-specific logits from a FROZEN normalized state representation.

    A separate fitted head is required for each declared prediction horizon.
    Its probability calibration and eligible risk sets belong to the evaluation
    pipeline; this class does not make an uncalibrated score a causal risk.
    """

    def __init__(self, state_dim: int, n_endpoints: int,
                 hidden_dim: int | Sequence[int] = (1024, 512),
                 dropout: float = 0.1):
        super().__init__()
        widths = [hidden_dim] if isinstance(hidden_dim, int) else list(hidden_dim)
        if min(int(state_dim), int(n_endpoints)) < 1 or not widths:
            raise ValueError("Risk-head dimensions must be positive")
        if any(int(width) != width or int(width) < 1 for width in widths):
            raise ValueError("Risk-head hidden dimensions must be positive integers")
        if not 0 <= float(dropout) < 1:
            raise ValueError("Risk-head dropout must be in [0,1)")
        self.n_endpoints = int(n_endpoints)
        self.hidden_dims = [int(width) for width in widths]
        layers: list[nn.Module] = []
        width = int(state_dim)
        for next_width in self.hidden_dims:
            layers.extend((nn.Linear(width, next_width), nn.GELU(),
                           nn.LayerNorm(next_width), nn.Dropout(float(dropout))))
            width = next_width
        layers.append(nn.Linear(width, self.n_endpoints))
        self.net = nn.Sequential(*layers)

    def forward(self, features: Tensor) -> Tensor:
        return self.net(features)
