"""Hierarchical molecular Transformer and differentiable person-to-person copying.

Outcome values never enter the encoder or query/key networks. During training,
all donors from the query's inner fold (and therefore its family) are masked.
Cached donor keys are refreshed each epoch, detached from the current gradient.
"""
from pathlib import Path
import copy
import time
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from common import log, dump
from attention import ContextStack, torch_mixture_scores


def device_for(request):
    if request == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if request == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable. Run --check-device or use --device cpu.")
    return request


def token_partition(x, features, count, seed, module_file=""):
    """Every retained assay belongs to exactly one token. No Y is used."""
    labels = np.full(len(features), -1, int)
    names = []
    if module_file:
        tab = pd.read_csv(module_file, sep="\t", dtype=str)
        if not {"feature", "module"} <= set(tab):
            raise ValueError("--module-file requires feature and module TSV columns")
        # Deterministic first alphabetical membership for overlapping pathways.
        lookup = {f: i for i, f in enumerate(features)}
        for name, sub in tab.sort_values(["module", "feature"]).groupby("module"):
            ix = [lookup[f] for f in sub.feature.unique() if f in lookup and labels[lookup[f]] < 0]
            if ix:
                labels[ix] = len(names); names.append(str(name))
    left = np.flatnonzero(labels < 0)
    if len(left):
        k = min(max(1, count-len(names)), len(left))
        if k == 1:
            cluster = np.zeros(len(left), int)
        else:
            d = min(24, len(x)-1, len(left))
            load = PCA(d, svd_solver="randomized", random_state=seed).fit(x[:, left]).components_.T
            load /= np.maximum(np.linalg.norm(load, axis=1, keepdims=True), 1e-8)
            cluster = MiniBatchKMeans(k, random_state=seed, n_init=5, batch_size=512).fit_predict(load)
        for group in np.unique(cluster):
            labels[left[cluster == group]] = len(names)
            names.append(f"data_token_{len(names)+1:03d}")
    return labels, names


class RowEncoder(nn.Module):
    def __init__(self, membership, width=64, heads=4, layers=2, dropout=.1,
                 architecture="transformer", cross_mode="qkv"):
        super().__init__()
        membership = np.asarray(membership, int)
        self.n_features, self.n_tokens = len(membership), int(membership.max())+1
        self.width, self.heads, self.architecture = width, heads, architecture
        self.cross_mode = cross_mode
        self.register_buffer("membership", torch.tensor(membership, dtype=torch.long))
        counts = np.bincount(membership, minlength=self.n_tokens).astype("float32")
        self.register_buffer("counts", torch.tensor(counts))
        # Learned feature identities distinguish every assay within its module.
        self.value_embedding = nn.Parameter(torch.randn(self.n_features, width)*.05)
        self.missing_embedding = nn.Parameter(torch.randn(self.n_features, width)*.02)
        self.token_embedding = nn.Parameter(torch.randn(self.n_tokens, width)*.02)
        self.cls = nn.Parameter(torch.randn(1, 1, width)*.02)
        self.input_norm = nn.LayerNorm(width)
        if architecture in ["transformer", "uniform"]:
            self.context = ContextStack(width, heads, layers, dropout,
                                        "uniform" if architecture == "uniform" else "learned")
        else:
            self.context = nn.Sequential(nn.Linear(self.n_tokens*width, width*4), nn.GELU(),
                nn.Dropout(dropout), nn.Linear(width*4, width), nn.LayerNorm(width))
        self.project = nn.Sequential(nn.Linear(width, width), nn.GELU(), nn.Linear(width, width))
        self.direct = nn.Linear(width, 1)
        self.decoder_weight = nn.Parameter(torch.randn(self.n_features, width)*.05)
        self.decoder_bias = nn.Parameter(torch.zeros(self.n_features))
        self.gate = nn.Linear(width, heads)
        self.temperature = nn.Parameter(torch.full((heads,), -1. if cross_mode=="qkv" else -3.))
        self.query = nn.Linear(width, width, bias=False)
        self.key = nn.Linear(width, width, bias=False)
        # An asymmetric learned Q/K projection; sqrt(width) restores unit scale
        # after the normalized row embedding, before scaled dot product.
        with torch.no_grad():
            self.query.weight.copy_(torch.eye(width)+torch.randn(width,width)*.01)
            self.key.weight.copy_(torch.eye(width)+torch.randn(width,width)*.01)

    def forward(self, x, observed, decode=False, return_attention=False):
        # No feature value hidden by observed=False can enter the encoder.
        x = torch.where(observed, x, torch.zeros_like(x))
        token = x.new_zeros((len(x), self.n_tokens, self.width))
        token.index_add_(1, self.membership, x.unsqueeze(-1)*self.value_embedding)
        token.index_add_(1, self.membership, (~observed).unsqueeze(-1)*self.missing_embedding)
        token = token/self.counts.sqrt()[None, :, None] + self.token_embedding
        token = self.input_norm(token)
        maps = None
        if self.architecture in ["transformer", "uniform"]:
            context, maps = self.context(torch.cat([self.cls.expand(len(x), -1, -1), token], dim=1), return_attention)
            row, contextual = context[:, 0], context[:, 1:]
        else:
            row = self.context(token.flatten(1))
            contextual = token + row[:, None, :]
        z = F.normalize(self.project(row), dim=-1)
        reconstruction = None
        if decode:
            reconstruction = (contextual[:, self.membership]*self.decoder_weight).sum(-1)+self.decoder_bias
        result = (z, self.direct(row).squeeze(-1), reconstruction)
        return (*result, maps) if return_attention else result

    def head_scores(self, q, bank):
        h, d = self.heads, self.width//self.heads
        scale = self.width**.5
        qh = self.query(q*scale).reshape(-1,h,d)
        kh = self.key(bank*scale).reshape(-1,h,d)
        temperature = .02 + 1.98*torch.sigmoid(self.temperature)
        return torch.einsum("bhd,rhd->brh", qh, kh)/(d**.5*temperature)

    def scores(self, q, bank):
        if self.cross_mode == "qkv":
            return torch_mixture_scores(self.head_scores(q,bank),torch.softmax(self.gate(q),-1))
        h, d = self.heads, self.width//self.heads
        qh, bh = q.reshape(-1, h, d), bank.reshape(-1, h, d)
        distances = (qh.square().sum(-1)[:, None, :] + bh.square().sum(-1)[None, :, :]
                     - 2*torch.einsum("bhd,rhd->brh", qh, bh)).clamp_min(0)
        gate = torch.softmax(self.gate(q), -1)
        temperature = .02 + 1.98*torch.sigmoid(self.temperature)
        return -(distances*gate[:, None, :]/temperature).sum(-1)


def torch_copy(model, q, bank, by, bw, forbidden, k, prior, strength):
    if torch.any((~forbidden).sum(1) < k):
        raise ValueError("Insufficient independent neural donors")
    q, bank = q.float(), bank.float()
    if model.cross_mode == "qkv":
        logits = model.head_scores(q,bank).masked_fill(forbidden[:,:,None],-torch.inf)
        gate = torch.softmax(model.gate(q),-1)
        score = torch_mixture_scores(logits,gate)
        jj = torch.argsort(score,descending=True,stable=True)[:,:k]
        chosen = logits.gather(1,jj[:,:,None].expand(-1,-1,model.heads))
        # Label donors only; IPCW is a positive measure over keys, not part of Q/K.
        head_weights = torch.softmax(chosen+bw[jj].log()[:,:,None],dim=1)
        weights = (head_weights*gate[:,None,:]).sum(-1)
    else:
        score = model.scores(q,bank).masked_fill(forbidden,-torch.inf)
        jj = torch.argsort(score,descending=True,stable=True)[:,:k]
        weights = torch.softmax(score.gather(1,jj),dim=1)*bw[jj]
        weights = weights/weights.sum(1,keepdim=True)
    ess = weights.square().sum(1).reciprocal()
    p = (ess*(weights*by[jj]).sum(1)+strength*prior)/(ess+strength)
    return p, jj, weights


@torch.no_grad()
def encode(model, x, observed, device="cpu", batch_size=256, decode=False):
    model.to(device).eval()
    zs, ps, rs = [], [], []
    for start in range(0, len(x), batch_size):
        xx = torch.as_tensor(np.ascontiguousarray(x[start:start+batch_size]), dtype=torch.float32, device=device)
        mm = torch.as_tensor(np.ascontiguousarray(observed[start:start+batch_size]), dtype=torch.bool, device=device)
        z, p, r = model(xx, mm, decode)
        zs.append(z.cpu().numpy()); ps.append(torch.sigmoid(p).cpu().numpy())
        if decode:
            rs.append(r.cpu().numpy())
    return np.concatenate(zs), np.concatenate(ps), np.concatenate(rs) if decode else None


def train_encoder(x, observed, y, w, ids, groups, membership, a, output, architecture="transformer",
                  objective="retrieval", pretrain=True, cross_mode="qkv"):
    """Only build rows enter here; an inner fold is reserved for early stopping.

    The OOF quality score used elsewhere is independently cross-fitted. This
    neural model's training predictions are never relabelled as OOF evidence.
    """
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    dev = device_for(a.device)
    torch.set_num_threads(a.cores)
    torch.manual_seed(a.seed); np.random.seed(a.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(a.seed)
    torch.use_deterministic_algorithms(a.deterministic, warn_only=False)
    cv = (StratifiedGroupKFold(5, shuffle=True, random_state=a.seed+43) if groups is not None
          else StratifiedKFold(5, shuffle=True, random_state=a.seed+43))
    strata = np.where(w > 0, y, 2)
    tr, va = next(cv.split(x, strata, groups))
    folds = np.full(len(x), -1, int)
    cv2 = (StratifiedGroupKFold(a.folds, shuffle=True, random_state=a.seed+44) if groups is not None
           else StratifiedKFold(a.folds, shuffle=True, random_state=a.seed+44))
    for fold, (_, ix) in enumerate(cv2.split(x[tr], strata[tr], None if groups is None else groups[tr])):
        folds[tr[ix]] = fold
    pd.DataFrame({"eid": ids, "inner_role": np.where(folds < 0, "early_stop", "train"),
                  "donor_exclusion_fold": folds}).to_csv(output/"inner_split.csv", index=False)
    donor = tr[w[tr] > 0]
    if min(np.bincount(y[donor], minlength=2)) < 5:
        raise ValueError("Neural inner training needs >=5 known outcomes in each class")
    k = min(a.neural_k, min(np.sum(folds[donor] != f) for f in np.unique(folds[tr])))
    if k < 2:
        raise ValueError("Too few fold-disjoint donors")
    model = RowEncoder(membership, a.width, a.heads, a.layers, a.dropout, architecture, cross_mode).to(dev)
    prior = float(np.average(y[donor], weights=w[donor]))
    with torch.no_grad():
        model.direct.bias.fill_(np.log(prior/(1-prior)))
    optimizer = torch.optim.AdamW(model.parameters(), lr=a.learning_rate, weight_decay=a.weight_decay)
    amp = a.amp and dev == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    rng = np.random.default_rng(a.seed+47)
    logs, completed, best, bad, best_state = [], {}, np.inf, 0, None
    checkpoint = output/"checkpoint.pt"
    ssl_state = None
    if a.resume and checkpoint.exists():
        saved = torch.load(checkpoint, map_location="cpu", weights_only=False)
        model.load_state_dict(saved["state"]); optimizer.load_state_dict(saved["optimizer"])
        for state in optimizer.state.values():
            for key, value in state.items():
                if torch.is_tensor(value):
                    state[key] = value.to(dev)
        scaler.load_state_dict(saved["scaler"])
        logs, completed, best, bad = saved["logs"], saved["completed"], saved["best"], saved["bad"]
        best_state, ssl_state = saved["best_state"], saved["ssl_state"]
        rng.bit_generator.state = saved["numpy_rng"]
        torch.set_rng_state(saved["torch_rng"])
        if dev == "cuda" and saved["cuda_rng"] is not None:
            torch.cuda.set_rng_state_all(saved["cuda_rng"])
        log("RESUME", output.name, str(completed))
    val_mask = observed[va] & (np.random.default_rng(a.seed+49).random(observed[va].shape)<a.mask_fraction)
    by = torch.tensor(y[donor], dtype=torch.float32, device=dev)
    bw = torch.tensor(w[donor], dtype=torch.float32, device=dev)
    bank_fold = torch.tensor(folds[donor], device=dev)
    start_time = time.monotonic()
    for phase, epochs in [("pretrain", a.pretrain_epochs if pretrain else 0), ("joint", a.epochs)]:
        if phase == "joint" and completed.get("joint", 0) == 0:
            best, bad, best_state = np.inf, 0, None
        if completed.get(phase, 0) >= epochs or (phase == "joint" and bad >= a.patience):
            continue
        for epoch in range(completed.get(phase, 0), epochs):
            epoch_started = time.monotonic()
            if dev == "cuda": torch.cuda.reset_peak_memory_stats()
            if phase == "joint":
                bz = torch.tensor(encode(model, x[donor], observed[donor], dev, a.batch_size)[0], device=dev)
            model.train()
            totals = np.zeros(4); nb = 0
            order = rng.permutation(tr)
            for begin in range(0, len(order), a.batch_size):
                ix = order[begin:begin+a.batch_size]
                xx = torch.tensor(x[ix], dtype=torch.float32, device=dev)
                mm = torch.tensor(observed[ix], dtype=torch.bool, device=dev)
                hidden = mm & (torch.rand(mm.shape, device=dev) < a.mask_fraction)
                optimizer.zero_grad(set_to_none=True)
                with torch.autocast(device_type=dev, dtype=torch.float16, enabled=amp):
                    z, direct, reconstruction = model(xx, mm & ~hidden, decode=True)
                    rec = (reconstruction.float()-xx).square()[hidden].mean() if hidden.any() else reconstruction.sum()*0
                    direct_loss = rec*0; retrieval_loss = rec*0
                    total = rec if phase == "pretrain" else a.reconstruction_weight*rec
                    if phase == "joint":
                        yy = torch.tensor(y[ix], dtype=torch.float32, device=dev)
                        ww = torch.tensor(w[ix], dtype=torch.float32, device=dev)
                        direct_loss = (F.binary_cross_entropy_with_logits(direct.float(), yy, reduction="none")*ww).sum()/ww.sum().clamp_min(1)
                        if objective == "retrieval":
                            forbidden = torch.tensor(folds[ix], device=dev)[:, None] == bank_fold[None, :]
                            # COPY is already a probability, not a neural logit.
                            # Keep matching and its BCE in float32: CUDA autocast
                            # rejects probability BCE and half precision rounds
                            # probabilities near 1 before the logarithm.
                            with torch.autocast(device_type=dev, enabled=False):
                                borrowed, _, _ = torch_copy(model, z.float(), bz.float(), by, bw,
                                                           forbidden, k, prior, a.train_prior)
                                retrieval_loss = (F.binary_cross_entropy(borrowed.clamp(1e-5, 1-1e-5), yy, reduction="none")*ww).sum()/ww.sum().clamp_min(1)
                            total = total+retrieval_loss+a.direct_weight*direct_loss
                        else:
                            total = total+direct_loss
                if not torch.isfinite(total):
                    raise FloatingPointError("Nonfinite neural loss; inspect input/learning rate")
                scaler.scale(total).backward()
                scaler.unscale_(optimizer)
                nn.utils.clip_grad_norm_(model.parameters(), a.gradient_clip)
                scaler.step(optimizer); scaler.update()
                totals += [float(v.detach()) for v in [total, rec, direct_loss, retrieval_loss]]; nb += 1
                if time.monotonic()-start_time > 55:
                    log("RUNNING", output.name, f"{phase} epoch={epoch+1}, batch={nb}")
                    start_time = time.monotonic()
            model.eval()
            with torch.no_grad():
                zv, pv, rv = encode(model, x[va], observed[va] & ~val_mask, dev, a.batch_size, True)
                val_rec = float(np.mean((rv[val_mask]-x[va][val_mask])**2))
                if phase == "pretrain":
                    score = val_rec; val_ll = np.nan
                else:
                    if objective == "retrieval":
                        bz = torch.tensor(encode(model, x[donor], observed[donor], dev, a.batch_size)[0], device=dev)
                        vals = []
                        for begin in range(0, len(va), a.batch_size):
                            q = torch.tensor(zv[begin:begin+a.batch_size], device=dev)
                            pred, _, _ = torch_copy(model, q, bz, by, bw,
                                torch.zeros((len(q), len(donor)), dtype=torch.bool, device=dev), k, prior, a.train_prior)
                            vals.append(pred.cpu().numpy())
                        pv = np.concatenate(vals)
                    pv = np.clip(pv, 1e-6, 1-1e-6)
                    val_ll = float(np.average(-(y[va]*np.log(pv)+(1-y[va])*np.log1p(-pv)), weights=w[va]))
                    score = val_ll
            improved = score < best-1e-6
            if improved:
                best, bad = score, 0
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
            completed[phase] = epoch+1
            row = dict(phase=phase, epoch=epoch+1, train_loss=totals[0]/nb, train_masked_mse=totals[1]/nb,
                train_direct_loss=totals[2]/nb, train_retrieval_loss=totals[3]/nb,
                early_stop_logloss=val_ll, early_stop_masked_mse=val_rec, best=improved,
                seconds=time.monotonic()-epoch_started,
                cuda_peak_allocated_mb=torch.cuda.max_memory_allocated()/2**20 if dev=="cuda" else None)
            logs.append(row); pd.DataFrame(logs).to_csv(output/"learning_curve.csv", index=False)
            if phase == "pretrain" and epoch+1 == epochs:
                # Start supervised training from the best masked-reconstruction epoch.
                model.load_state_dict(best_state)
                ssl_state = copy.deepcopy(best_state)
            saved = dict(state={k:v.detach().cpu() for k,v in model.state_dict().items()},
                optimizer=optimizer.state_dict(), scaler=scaler.state_dict(), logs=logs, completed=completed,
                best=best, bad=bad, best_state=best_state, ssl_state=ssl_state,
                numpy_rng=rng.bit_generator.state, torch_rng=torch.get_rng_state(),
                cuda_rng=torch.cuda.get_rng_state_all() if dev == "cuda" else None)
            tmp = checkpoint.with_suffix(".tmp"); torch.save(saved, tmp); tmp.replace(checkpoint)
            logloss_text = "N/A (reconstruction only)" if phase == "pretrain" else f"{val_ll:.4f}"
            log("EPOCH", output.name, f"{phase} {epoch+1}/{epochs}; masked_MSE={val_rec:.4f}; logloss={logloss_text}; seconds={row['seconds']:.1f}")
            if phase == "joint" and bad >= a.patience:
                break
    if best_state is None:
        raise RuntimeError("No neural checkpoint selected")
    model.load_state_dict(best_state); model.cpu().eval()
    ssl = copy.deepcopy(model) if ssl_state is not None else None
    if ssl is not None:
        ssl.load_state_dict(ssl_state)
    dump(output/"training_summary.json", dict(device=dev, architecture=architecture, objective=objective,
        parameters=sum(p.numel() for p in model.parameters()), build=len(x), gradient_rows=len(tr),
        early_stop_rows=len(va), known_donors=len(donor), donor_fold_exclusion=True,
        bank_refresh="each_epoch; detached row embeddings; Q/K projections remain differentiable", cross_mode=cross_mode, completed_epochs=completed,
        best_early_stop_score=best, torch_version=torch.__version__, mixed_precision=amp))
    torch.save(model.state_dict(), output/"selected_weights.pt")
    return model, ssl
