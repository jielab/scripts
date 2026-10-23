"""CSx block-network predictor for the FINAL ABM workflow."""
from pathlib import Path
import math
import time

import numpy as np
import pandas as pd
import torch
from torch import nn
from scipy.special import expit

def available_memory():
    """Conservative host/cgroup headroom; zero disables RAM staging on unknown hosts."""
    try:
        fields=dict(line.split(":",1) for line in Path("/proc/meminfo").read_text().splitlines())
        available=int(fields["MemAvailable"].split()[0])*1024
        # Resolve this process's cgroup, including nested container limits.
        cgroup=next((line.split(":",2)[2] for line in Path("/proc/self/cgroup").read_text().splitlines()
                     if line.startswith("0::")), "")
        root=Path("/sys/fs/cgroup")
        current=root/cgroup.lstrip("/")
        while current.is_relative_to(root):
            limit=current/"memory.max"
            if limit.is_file() and limit.read_text().strip()!="max":
                available=min(available,max(0,int(limit.read_text())-int((current/"memory.current").read_text())))
            if current==root:break
            current=current.parent
        return available
    except (OSError,ValueError,KeyError):
        return 0

def stage_block_matrix(a,cfg):
    """Sequentially copy a single trait's matrix; preserve raw values and row order."""
    limit=float(cfg.get("block_ram_cache_gib",16))*2**30
    if not np.isfinite(limit) or limit<0:raise ValueError("block_ram_cache_gib must be finite and nonnegative")
    budget=min(limit,available_memory()//4)
    if a.nbytes>budget or not limit:
        print(f"[GRID ABM] REUSE block matrix via mmap: {a.nbytes/2**30:.2f} GiB; RAM budget={budget/2**30:.2f} GiB",flush=True)
        return a
    print(f"[GRID ABM] START block matrix -> RAM: {a.nbytes/2**30:.2f} GiB",flush=True)
    try:
        cached=np.empty(a.shape,dtype=a.dtype)
        batch=max(1,(128*1024**2)//max(1,a[0:1].nbytes))
        last=time.monotonic()
        for start in range(0,len(a),batch):
            end=min(start+batch,len(a));cached[start:end]=a[start:end]
            now=time.monotonic()
            if now-last>=60:
                print(f"[GRID ABM] PROGRESS block matrix -> RAM: {end:,}/{len(a):,}",flush=True);last=now
    except MemoryError:
        print("[GRID ABM] REUSE block matrix via mmap: RAM allocation unavailable",flush=True)
        return a
    cached.setflags(write=False)
    print("[GRID ABM] DONE block matrix -> RAM; reuse for scaling, all epochs and predictions",flush=True)
    return cached

def feature_scaler(a,rr,batch=512):
    n=0;s=np.zeros(a.shape[1:],np.float64);ss=np.zeros(a.shape[1:],np.float64)
    for st in range(0,len(rr),batch):
        x=np.asarray(a[rr[st:st+batch]],np.float32)
        if not np.isfinite(x).all():raise ValueError("Nonfinite block score")
        n+=len(x);s+=x.sum(0,dtype=np.float64);ss+=np.square(x,dtype=np.float64).sum(0,dtype=np.float64)
    mu=s/n;var=np.maximum(ss/n-mu*mu,1e-12)
    return mu.astype(np.float32),np.sqrt(var).astype(np.float32)

class BlockTransformer(nn.Module):
    """LD-block PRS tokens -> latent genetic embedding -> phenotype head.

    Declared side features enter the prediction head. ABM support uses separate fixed genetic geometry.
    """
    def __init__(self,tokens,channels,side_dim,cfg):
        super().__init__();d=int(cfg["d_model"]);self.channels=channels
        self.value=nn.Linear(channels,d,bias=False)
        self.token=nn.Embedding(len(tokens),d)
        self.chr=nn.Embedding(22,d)
        popmap={"AFR":0,"EAS":1,"EUR":2,"SAS":3,"ALL":4}
        self.pop=nn.Embedding(5,d)
        self.register_buffer("token_idx",torch.arange(len(tokens),dtype=torch.long))
        self.register_buffer("chr_idx",torch.tensor(tokens.CHR.to_numpy(int)-1,dtype=torch.long))
        self.register_buffer("pop_idx",torch.tensor([popmap[x] for x in tokens.source_population],dtype=torch.long))
        cont=np.c_[tokens.midpoint.to_numpy(float),np.log1p(tokens.length.to_numpy(float))]
        cont=(cont-cont.mean(0))/np.maximum(cont.std(0),1)
        self.register_buffer("cont",torch.tensor(cont,dtype=torch.float32))
        self.meta=nn.Sequential(nn.Linear(2,d),nn.GELU(),nn.Linear(d,d))
        self.norm=nn.LayerNorm(d)
        self.latent=nn.Parameter(torch.empty(int(cfg["latent_tokens"]),d));nn.init.normal_(self.latent,std=.02)
        self.cross=nn.MultiheadAttention(d,int(cfg["heads"]),dropout=float(cfg["dropout"]),batch_first=True)
        self.blocks=nn.ModuleList([nn.TransformerEncoderLayer(d,int(cfg["heads"]),4*d,float(cfg["dropout"]),
                            activation="gelu",batch_first=True,norm_first=True) for _ in range(int(cfg["layers"]))])
        self.out=nn.LayerNorm(d)
        self.head=nn.Sequential(nn.Linear(d+side_dim,d),nn.GELU(),nn.Dropout(float(cfg["dropout"])),nn.Linear(d,1))
    def forward(self,x,side):
        z=self.value(x)+self.token(self.token_idx)[None]+self.chr(self.chr_idx)[None]+self.pop(self.pop_idx)[None]+self.meta(self.cont)[None]
        z=self.norm(z);q=self.latent[None].expand(len(x),-1,-1)
        h=q+self.cross(q,z,z,need_weights=False)[0]
        for b in self.blocks:h=b(h)
        emb=self.out(h.mean(1))
        return self.head(torch.cat([emb,side],1)).squeeze(1),emb

class BlockLinear(nn.Module):
    def __init__(self,tokens,channels,side_dim):super().__init__();self.head=nn.Linear(tokens*channels+side_dim,1)
    def forward(self,x,side):return self.head(torch.cat([x.flatten(1),side],1)).squeeze(1)

def batch_data(a,rr,ix,mu,sd,side,device):
    x=(np.asarray(a[rr[ix]],np.float32)-mu)/sd
    return torch.as_tensor(x,device=device),torch.as_tensor(side[ix],device=device)

def batches(ix,n,seed=None):
    z=np.asarray(ix,np.int64).copy()
    if seed is not None:np.random.default_rng(seed).shuffle(z)
    for st in range(0,len(z),n):yield z[st:st+n]

def row_loss(y,p,binary):
    if binary:
        p=np.clip(p,1e-7,1-1e-7);return -(y*np.log(p)+(1-y)*np.log1p(-p))
    return np.square(y-p)

def predict_model(model,kind,a,rr,side,ix,mu,sd,cfg,device,ymean,ysd,embedding=False):
    model.eval();pp=[];ee=[];seen=0;last_progress=time.monotonic()
    with torch.inference_mode():
        for b in batches(ix,int(cfg["batch_size"])):
            x,s=batch_data(a,rr,b,mu,sd,side,device)
            with torch.autocast(device.type,dtype=torch.bfloat16,enabled=device.type=="cuda" and cfg.get("amp")=="bf16"):
                if kind=="transformer":raw,e=model(x,s)
                else:raw=model(x,s);e=None
            r=raw.float().cpu().numpy()
            if not cfg["_binary"]:r=r*ysd+ymean
            pp.append(r)
            if embedding:ee.append(e.float().cpu().numpy())
            seen+=len(b);now=time.monotonic()
            if now-last_progress>=60:
                print(f"[GRID ABM] PROGRESS predict {kind}: {seen:,}/{len(ix):,}",flush=True);last_progress=now
    return np.concatenate(pp),np.concatenate(ee) if embedding else None

def train(kind,a,rr,side,y,roles,tokens,cfg,device,mu,sd):
    bi=np.flatnonzero(roles=="build");tu=np.flatnonzero(roles=="tune_model")
    model=BlockTransformer(tokens,a.shape[2],side.shape[1],cfg).to(device) if kind=="transformer" else BlockLinear(a.shape[1],a.shape[2],side.shape[1]).to(device)
    opt=torch.optim.AdamW(model.parameters(),lr=float(cfg["learning_rate"]),weight_decay=float(cfg["weight_decay"]))
    binary=cfg["_binary"];ymean=None;ysd=None
    if binary:
        cases=max(np.sum(y[bi]==1),1);controls=max(np.sum(y[bi]==0),1)
        if cfg.get("initialize_prevalence"):
            with torch.no_grad():
                model.head[-1].bias.fill_(math.log(cases/controls))
        if cfg.get("binary_objective")=="brier":
            criterion=lambda logits,target: torch.mean((torch.sigmoid(logits)-target)**2)
        else:
            criterion=nn.BCEWithLogitsLoss(pos_weight=torch.tensor(controls/cases,dtype=torch.float32,device=device))
    else:
        ymean=float(np.mean(y[bi]));ysd=max(float(np.std(y[bi])),1e-8);criterion=nn.MSELoss()
    print(f"[GRID ABM] START {kind}: device={next(model.parameters()).device}, train={len(bi):,}, tune={len(tu):,}, epochs={cfg['epochs']}",flush=True)
    best=None;hist=[];pat=0
    for epoch in range(1,int(cfg["epochs"])+1):
        model.train();ls=[];seen=0;last_progress=time.monotonic()
        for ix in batches(bi,int(cfg["batch_size"]),int(cfg["seed"])+epoch+(1000 if kind=="linear" else 0)):
            x,s=batch_data(a,rr,ix,mu,sd,side,device);yy=torch.as_tensor(y[ix],dtype=torch.float32,device=device)
            if not binary:yy=(yy-ymean)/ysd
            opt.zero_grad(set_to_none=True)
            with torch.autocast(device.type,dtype=torch.bfloat16,enabled=device.type=="cuda" and cfg.get("amp")=="bf16"):
                raw=model(x,s)[0] if kind=="transformer" else model(x,s);loss=criterion(raw.float(),yy)
            if not torch.isfinite(loss):raise RuntimeError("Nonfinite blocknet training loss")
            loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),float(cfg["gradient_clip"]));opt.step();ls.append(loss.detach())
            seen+=len(ix)
            now=time.monotonic()
            if len(ls)==1 or seen==len(bi) or now-last_progress>=60:
                print(f"[GRID ABM] PROGRESS {kind} epoch={epoch}/{cfg['epochs']} train={seen:,}/{len(bi):,} device={x.device}",flush=True)
                last_progress=now
        pv,_=predict_model(model,kind,a,rr,side,tu,mu,sd,cfg,device,ymean,ysd)
        tl=(float(np.mean((y[tu]-expit(pv))**2)) if binary and cfg.get("binary_objective")=="brier"
            else float(np.mean(row_loss(y[tu],expit(pv) if binary else pv,binary))))
        hist.append(dict(epoch=epoch,train_loss=float(np.mean(torch.stack(ls).cpu().numpy().astype(np.float64))),tune_loss=tl))
        print(f"[GRID ABM] {kind} epoch={epoch} tune_loss={tl:.8g}",flush=True)
        if best is None or tl<best[0]-float(cfg["min_delta"]):
            best=(tl,{k:v.detach().cpu().clone() for k,v in model.state_dict().items()},epoch);pat=0
        else:
            pat+=1
            if pat>=int(cfg["patience"]):break
    model.load_state_dict(best[1]);return model,pd.DataFrame(hist),dict(best_epoch=best[2],ymean=ymean,ysd=ysd)
