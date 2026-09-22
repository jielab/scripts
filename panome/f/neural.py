"""Mini-batch neural prediction with an exact piecewise-exponential likelihood.

Censored participants contribute their observed exposure time. This avoids the
incorrect practice of computing a Cox partial likelihood on incomplete batch risk sets.
"""
import copy
import hashlib
import numpy as np
import pandas as pd
import torch
from torch import nn
from common import log
from representation import device_for, initialize

def exposure_matrix(times, cuts):
    times = np.asarray(times, float)
    left = np.r_[0, cuts]
    width = np.diff(np.r_[left, np.inf])
    return np.minimum(np.maximum(times[:, None] - left, 0), width).astype("float32")

def random_neighbors(atlas, ids, groups, seed):
    result = np.empty((len(ids), atlas.k), dtype=int)
    for row, person in enumerate(ids):
        hashed = int.from_bytes(hashlib.sha256(f"{seed}:{person}".encode()).digest()[:8], "little")
        rng = np.random.default_rng(hashed)
        selected = []
        while len(selected) < atlas.k:
            draw = rng.choice(len(atlas.train_ids), min(len(atlas.train_ids), atlas.k*3), replace=False)
            valid = atlas.train_ids[draw] != str(person)
            if groups is not None and atlas.train_groups is not None:
                valid &= atlas.train_groups[draw] != str(groups[row])
            selected = list(dict.fromkeys(selected + draw[valid].tolist()))
        result[row] = selected[:atlas.k]
    return result

class RiskNetwork(nn.Module):
    def __init__(self, n_input, n_clinical, mode="mlp", dim=32, heads=4):
        super().__init__()
        self.config = dict(n_input=n_input, n_clinical=n_clinical, mode=mode,
                           dim=dim, heads=heads)
        self.mode = mode
        if mode == "row":
            self.query = nn.Linear(n_input, dim)
            self.memory = nn.Linear(n_input, dim)
            self.attention = nn.MultiheadAttention(dim, heads, dropout=.1, batch_first=True)
            size = n_input + dim + n_clinical
        else:
            size = n_input + n_clinical
        self.risk = nn.Sequential(nn.Linear(size, 64), nn.GELU(),
                                  nn.Dropout(.1), nn.Linear(64, 1))

    def forward(self, x, clinical, references=None, return_weights=False):
        attention = None
        if self.mode == "row":
            query = self.query(x).unsqueeze(1)
            memory = self.memory(references)
            context, attention = self.attention(query, memory, memory,
                                                need_weights=return_weights)
            x = torch.cat([x, context.squeeze(1)], 1)
        result = self.risk(torch.cat([clinical, x], 1)).squeeze(1)
        return (result, attention) if return_weights else result

class NeuralModel:
    def __init__(self, mode, outcome_type, dim=32, heads=4):
        self.mode, self.outcome_type, self.dim, self.heads = mode, outcome_type, dim, heads

    def fit(self, x, clinical, p, neighbors, reference_bank, a, out, name):
        initialize(a)
        device = device_for(a)
        self.network = RiskNetwork(x.shape[1], clinical.shape[1], self.mode, self.dim, self.heads).to(device)
        train = np.flatnonzero(p.split.to_numpy() == "train")
        valid = np.flatnonzero(p.split.to_numpy() == "validation")
        params = list(self.network.parameters())
        if self.outcome_type == "survival":
            times = p.time.to_numpy(float)
            event = p.event.to_numpy(float)
            event_times = times[train][event[train].astype(bool)]
            self.cuts = np.unique(np.quantile(event_times, np.linspace(0, 1, a.time_bins+1)[1:-1]))
            self.cuts = self.cuts[self.cuts > 0]
            exposure = exposure_matrix(times, self.cuts)
            event_bin = np.searchsorted(self.cuts, times, side="right")
            counts = np.bincount(event_bin[train], weights=event[train], minlength=len(self.cuts)+1)
            rates = (counts+.5)/(exposure[train].sum(0)+1)
            base = nn.Parameter(torch.as_tensor(np.log(rates), dtype=torch.float32, device=device))
            params.append(base)
        else:
            y = p.target.to_numpy(float)
            self.target_mean, self.target_scale = float(y[train].mean()), float(y[train].std())
            if self.target_scale < 1e-8:
                raise ValueError("Constant quantitative training target")
            target = (y-self.target_mean)/self.target_scale
            self.cuts = np.array([])
            base = None
        optimizer = torch.optim.AdamW(params, lr=a.neural_lr, weight_decay=a.neural_weight_decay)
        rng = np.random.default_rng(a.seed)
        best, wait, checkpoint = np.inf, 0, None
        history = []

        def loss_for(ix):
            xx = torch.as_tensor(np.array(x[ix]), device=device)
            cc = torch.as_tensor(np.array(clinical[ix]), device=device)
            refs = None
            if self.mode == "row":
                refs = torch.as_tensor(reference_bank[neighbors[ix]], device=device)
            score = self.network(xx, cc, refs)
            if self.outcome_type == "survival":
                exp_time = torch.as_tensor(exposure[ix], device=device)
                log_rates = (score[:, None] + base[None, :]).clamp(-25, 20)
                cumulative = (torch.exp(log_rates) * exp_time).sum(1)
                hazard_event = log_rates[torch.arange(len(ix), device=device),
                                        torch.as_tensor(event_bin[ix], device=device)]
                loss = (cumulative - torch.as_tensor(event[ix], dtype=torch.float32, device=device)*hazard_event).mean()
            else:
                loss = (score-torch.as_tensor(target[ix], dtype=torch.float32, device=device)).square().mean()
            return loss

        for epoch in range(a.neural_epochs):
            self.network.train()
            order = rng.permutation(train)
            total = 0
            for start in range(0, len(order), a.batch_size):
                ix = order[start:start+a.batch_size]
                loss = loss_for(ix)
                if not torch.isfinite(loss):
                    raise ValueError(f"{name}: nonfinite prediction loss")
                optimizer.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(params, 5)
                optimizer.step()
                total += loss.item()*len(ix)
            self.network.eval()
            vn = 0
            with torch.no_grad():
                for start in range(0, len(valid), a.batch_size):
                    ix = valid[start:start+a.batch_size]
                    vn += loss_for(ix).item()*len(ix)
            value = vn/len(valid)
            history.append(dict(epoch=epoch+1, train_loss=total/len(train), validation_loss=value))
            if np.isfinite(value) and value < best-1e-5:
                best, wait = value, 0
                checkpoint = (copy.deepcopy(self.network.state_dict()),
                              None if base is None else base.detach().cpu().numpy().copy())
            else:
                wait += 1
            if epoch == 0 or (epoch+1) % a.log_every == 0 or wait >= a.patience:
                log("FIT", name, f"epoch={epoch+1}, validation_loss={value:.5f}")
            if wait >= a.patience:
                break
        if checkpoint is None:
            raise ValueError(f"{name}: no finite checkpoint")
        self.network.load_state_dict(checkpoint[0])
        self.network.cpu().eval()
        self.base_log_rates = checkpoint[1]
        self.best_validation_loss = best
        pd.DataFrame(history).to_csv(out / f"{name}_history.csv", index=False)
        torch.save(dict(config=self.network.config, outcome_type=self.outcome_type,
            state_dict={k:v.cpu() for k,v in checkpoint[0].items()},
            cuts=self.cuts.tolist(),
            base_log_rates=None if checkpoint[1] is None else checkpoint[1].tolist(),
            target_mean=getattr(self, "target_mean", None),
            target_scale=getattr(self, "target_scale", None)), out / f"{name}.pt")
        return self

    def predict(self, x, clinical, neighbors=None, reference_bank=None, batch=512, return_attention=False):
        self.network.eval()
        predictions, attentions = [], []
        with torch.no_grad():
            for start in range(0, len(x), batch):
                stop = start+batch
                xx = torch.as_tensor(np.array(x[start:stop]), dtype=torch.float32)
                cc = torch.as_tensor(np.array(clinical[start:stop]), dtype=torch.float32)
                refs = None if self.mode != "row" else torch.as_tensor(
                    reference_bank[neighbors[start:stop]], dtype=torch.float32)
                score, weights = self.network(xx, cc, refs, return_weights=True)
                predictions.append(score.numpy())
                if self.mode == "row" and return_attention:
                    attentions.append(weights.squeeze(1).numpy())
        prediction = np.concatenate(predictions)
        if self.outcome_type == "quantitative":
            prediction = prediction*self.target_scale + self.target_mean
        return (prediction, np.concatenate(attentions)) if attentions else prediction

    def survival(self, score, horizon):
        rates = np.exp(np.clip(np.asarray(score)[:, None] + self.base_log_rates[None, :], -25, 20))
        cumulative = rates @ exposure_matrix([horizon], self.cuts)[0]
        return np.exp(-cumulative)

    @classmethod
    def load(cls, path):
        ck = torch.load(path, map_location="cpu", weights_only=True)
        model = cls(ck["config"]["mode"], ck["outcome_type"],
                    ck["config"]["dim"], ck["config"]["heads"])
        model.network = RiskNetwork(**ck["config"])
        model.network.load_state_dict(ck["state_dict"])
        model.network.eval()
        model.cuts = np.asarray(ck["cuts"])
        model.base_log_rates = None if ck["base_log_rates"] is None else np.asarray(ck["base_log_rates"])
        model.target_mean, model.target_scale = ck["target_mean"], ck["target_scale"]
        return model
