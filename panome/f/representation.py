"""Denoising autoencoder with observed-entry loss and validation early stopping."""
import copy
import numpy as np
import pandas as pd
import torch
from torch import nn

class Autoencoder(nn.Module):
    def __init__(self,p,latent,hidden):
        super().__init__()
        self.encoder=nn.Sequential(nn.Linear(p,hidden),nn.GELU(),nn.Linear(hidden,64),nn.GELU(),nn.Linear(64,latent))
        self.decoder=nn.Sequential(nn.Linear(latent,64),nn.GELU(),nn.Linear(64,hidden),nn.GELU(),nn.Linear(hidden,p))
    def forward(self,x): return self.decoder(self.encoder(x))


def encode(model,x,device,batch=512):
    model.eval()
    with torch.no_grad():
        return np.concatenate([model.encoder(torch.as_tensor(np.array(x[i:i+batch]),device=device)).cpu().numpy()
                               for i in range(0,len(x),batch)])


def fit_ae(x,mask,part,a,out):
    torch.manual_seed(a.seed)
    torch.set_num_threads(a.cores)
    if a.device=='cuda' and not torch.cuda.is_available(): raise ValueError('CUDA requested but unavailable')
    device='cuda' if a.device=='auto' and torch.cuda.is_available() else ('cpu' if a.device=='auto' else a.device)
    torch.use_deterministic_algorithms(True,warn_only=True)
    model=Autoencoder(x.shape[1],a.latent,a.hidden).to(device)
    optimizer=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
    tr=np.flatnonzero(part=='train'); va=np.flatnonzero(part=='validation')
    rng=np.random.default_rng(a.seed)
    best=np.inf; wait=0; history=[]
    for epoch in range(a.epochs):
        model.train(); numerator=denominator=0.
        for ix in np.array_split(rng.permutation(tr),max(1,int(np.ceil(len(tr)/a.batch_size)))):
            target=torch.as_tensor(x[ix],device=device)
            observed=torch.as_tensor(mask[ix],device=device)
            corrupted=target.clone()
            corrupted[torch.rand_like(corrupted)<a.corruption]=0
            loss=((model(corrupted)-target).square()*observed).sum()/observed.sum().clamp(min=1)
            optimizer.zero_grad();loss.backward();optimizer.step()
            numerator+=loss.item()*observed.sum().item();denominator+=observed.sum().item()
        model.eval(); vn=vd=0.
        with torch.no_grad():
            for start in range(0,len(va),a.batch_size):
                ix=va[start:start+a.batch_size]
                target=torch.as_tensor(x[ix],device=device);observed=torch.as_tensor(mask[ix],device=device)
                vn+=((model(target)-target).square()*observed).sum().item();vd+=observed.sum().item()
        value=vn/max(vd,1)
        history.append(dict(epoch=epoch+1,train_mse=numerator/denominator,validation_mse=value))
        if value<best-1e-5:best=value;state=copy.deepcopy(model.state_dict());wait=0
        else:wait+=1
        print(f'AE epoch {epoch+1}: validation MSE={value:.5f}',flush=True)
        if wait>=a.patience:break
    model.load_state_dict(state)
    torch.save(dict(state_dict={k:v.cpu() for k,v in state.items()},p=x.shape[1],latent=a.latent,hidden=a.hidden),out/'autoencoder.pt')
    pd.DataFrame(history).to_csv(out/'ae_history.csv',index=False)
    return encode(model,x,device),dict(device=device,best_validation_mse=best,epochs=len(history))
