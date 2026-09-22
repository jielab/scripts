"""Unsupervised encoders: endpoint-blind, frozen before prediction."""
import copy
import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.cluster import MiniBatchKMeans
from common import log

def device_for(a):
    if a.device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA requested but unavailable")
    return "cuda" if a.device == "auto" and torch.cuda.is_available() else ("cpu" if a.device == "auto" else a.device)

def initialize(a):
    torch.manual_seed(a.seed)
    np.random.seed(a.seed)
    torch.set_num_threads(a.cores)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(a.seed)
    torch.use_deterministic_algorithms(True, warn_only=True)

class Autoencoder(nn.Module):
    def __init__(self, p, latent=20, hidden=512):
        super().__init__()
        self.config = dict(p=p, latent=latent, hidden=hidden)
        self.encoder = nn.Sequential(nn.Linear(p * 2, hidden), nn.GELU(),
            nn.Linear(hidden, 64), nn.GELU(), nn.Linear(64, latent))
        self.decoder = nn.Sequential(nn.Linear(latent, 64), nn.GELU(),
            nn.Linear(64, hidden), nn.GELU(), nn.Linear(hidden, p))

    def encode(self, x, mask):
        return self.encoder(torch.cat([x * mask, mask], dim=1))

    def forward(self, x, mask):
        return self.decoder(self.encode(x, mask))

class MolecularTransformer(nn.Module):
    """Train-defined feature groups produce module tokens, not named pathways."""
    def __init__(self, p, groups, dim=32, heads=4, layers=2):
        super().__init__()
        self.config = dict(p=p, groups=groups, dim=dim, heads=heads, layers=layers)
        self.projections = nn.ModuleList()
        for j, indices in enumerate(groups):
            self.register_buffer(f"indices_{j}", torch.tensor(indices, dtype=torch.long))
            self.projections.append(nn.Linear(len(indices) * 2, dim))
        self.token_identity = nn.Parameter(torch.randn(len(groups), dim) * .02)
        layer = nn.TransformerEncoderLayer(dim, heads, dim * 2, dropout=.1,
                    activation="gelu", batch_first=True, norm_first=True)
        self.transformer = nn.TransformerEncoder(layer, layers, enable_nested_tensor=False)
        self.normalization = nn.LayerNorm(dim)
        self.decoder = nn.Linear(dim, p)

    def encode(self, x, mask):
        tokens = []
        for j, projection in enumerate(self.projections):
            ix = getattr(self, f"indices_{j}")
            tokens.append(projection(torch.cat([x[:, ix] * mask[:, ix], mask[:, ix]], 1)))
        tokens = torch.stack(tokens, 1) + self.token_identity
        return self.normalization(self.transformer(tokens).mean(1))

    def forward(self, x, mask):
        return self.decoder(self.encode(x, mask))

def groups_from_pca(pca, n_modules, seed):
    loadings = pca.components_.T.copy()
    loadings /= np.maximum(np.linalg.norm(loadings, axis=1, keepdims=True), 1e-8)
    count = min(n_modules, len(loadings))
    labels = MiniBatchKMeans(n_clusters=count, random_state=seed,
        n_init=5, batch_size=1024).fit_predict(loadings)
    return [np.flatnonzero(labels == j).tolist() for j in range(count)
            if np.any(labels == j)]

def encode(model, x, observed, device="cpu", batch=512):
    model = model.to(device).eval()
    values = []
    with torch.no_grad():
        for start in range(0, len(x), batch):
            xx = torch.as_tensor(np.array(x[start:start+batch]), device=device)
            mm = torch.as_tensor(np.array(observed[start:start+batch]), dtype=torch.float32, device=device)
            values.append(model.encode(xx, mm).cpu().numpy())
    return np.concatenate(values).astype("float32")

def fit_reconstruction(model, x, observed, part, a, out, name):
    initialize(a)
    device = device_for(a)
    model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=1e-4)
    train = np.flatnonzero(part == "train")
    valid = np.flatnonzero(part == "validation")
    rng = np.random.default_rng(a.seed)
    best, wait, saved = np.inf, 0, None
    history = []
    for epoch in range(a.epochs):
        model.train()
        numerator = denominator = 0
        order = rng.permutation(train)
        for start in range(0, len(order), a.batch_size):
            ix = order[start:start+a.batch_size]
            target = torch.as_tensor(np.array(x[ix]), device=device)
            mask = torch.as_tensor(np.array(observed[ix]), dtype=torch.float32, device=device)
            input_mask = mask * (torch.rand_like(mask) >= a.corruption)
            prediction = model(target, input_mask)
            loss = ((prediction - target).square() * mask).sum() / mask.sum().clamp(min=1)
            if not torch.isfinite(loss):
                raise ValueError(f"{name}: nonfinite reconstruction loss")
            opt.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5)
            opt.step()
            numerator += loss.item() * mask.sum().item()
            denominator += mask.sum().item()
        model.eval()
        vn = vd = 0
        with torch.no_grad():
            for start in range(0, len(valid), a.batch_size):
                ix = valid[start:start+a.batch_size]
                target = torch.as_tensor(np.array(x[ix]), device=device)
                mask = torch.as_tensor(np.array(observed[ix]), dtype=torch.float32, device=device)
                vn += (((model(target, mask) - target).square()) * mask).sum().item()
                vd += mask.sum().item()
        value = vn / max(vd, 1)
        history.append(dict(epoch=epoch+1, train_mse=numerator/max(denominator, 1),
                            validation_mse=value))
        if np.isfinite(value) and value < best - 1e-5:
            best, wait = value, 0
            saved = copy.deepcopy(model.state_dict())
        else:
            wait += 1
        if epoch == 0 or (epoch+1) % a.log_every == 0 or wait >= a.patience:
            log("FIT", name, f"epoch={epoch+1}, validation_MSE={value:.5f}")
        if wait >= a.patience:
            break
    if saved is None:
        raise ValueError(f"{name}: no finite validation checkpoint")
    model.load_state_dict(saved)
    torch.save(dict(kind=name, config=model.config,
                    state_dict={k:v.cpu() for k,v in saved.items()}), out / f"{name}.pt")
    pd.DataFrame(history).to_csv(out / f"{name}_history.csv", index=False)
    values = encode(model, x, observed, device, a.batch_size)
    return values, dict(best_validation_mse=best, epochs=len(history), device=device)

def load_encoder(path):
    checkpoint = torch.load(path, map_location="cpu", weights_only=True)
    cls = Autoencoder if checkpoint["kind"] == "ae" else MolecularTransformer
    model = cls(**checkpoint["config"])
    model.load_state_dict(checkpoint["state_dict"])
    return model.eval()
