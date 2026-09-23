"""Q/K/V attention with inspectable heads, and label-preserving cross attention.

Self attention has learned Q,K,V and output projections, residuals, pre-LN,
FFN and dropout. Cross attention uses learned Q,K but observed donor Y as V;
convex head mixing preserves a literal decomposition into real people's labels.
"""
import math
import numpy as np
from scipy.special import softmax, logsumexp
import torch
from torch import nn
from torch.nn import functional as F


class ContextBlock(nn.Module):
    def __init__(self, width, heads, dropout, mode="learned"):
        super().__init__()
        self.heads, self.width, self.mode = heads, width, mode
        self.norm1, self.norm2 = nn.LayerNorm(width), nn.LayerNorm(width)
        self.qkv = nn.Linear(width, 3*width)
        self.out = nn.Linear(width, width)
        self.ffn = nn.Sequential(nn.Linear(width, 4*width), nn.GELU(), nn.Dropout(dropout), nn.Linear(4*width, width))
        self.dropout = nn.Dropout(dropout)
        self.attention_dropout = dropout
        self.intervention = None  # evaluation-only intervention; never a fitted parameter

    def forward(self, x, capture=False):
        b, n, _ = x.shape
        q, k, v = self.qkv(self.norm1(x)).view(b,n,3,self.heads,self.width//self.heads).permute(2,0,3,1,4)
        mode = self.intervention or self.mode
        weights = None
        if mode == "learned" and not capture:
            context = F.scaled_dot_product_attention(q,k,v,dropout_p=self.attention_dropout if self.training else 0.)
        else:
            if mode == "uniform":
                weights = x.new_full((b,self.heads,n,n), 1/n)
            elif mode == "identity":
                weights = torch.eye(n,device=x.device,dtype=x.dtype)[None,None].expand(b,self.heads,-1,-1)
            else:
                weights = torch.softmax(q@k.transpose(-2,-1)/math.sqrt(q.shape[-1]),dim=-1)
            context = F.dropout(weights,p=self.attention_dropout,training=self.training)@v
        x = x+self.dropout(self.out(context.transpose(1,2).reshape(b,n,self.width)))
        x = x+self.dropout(self.ffn(self.norm2(x)))
        return x, weights


class ContextStack(nn.Module):
    def __init__(self,width,heads,layers,dropout,mode="learned"):
        super().__init__()
        self.layers = nn.ModuleList([ContextBlock(width,heads,dropout,mode) for _ in range(layers)])
        self.norm = nn.LayerNorm(width)

    def forward(self,x,capture=False):
        maps = []
        for layer in self.layers:
            x, weights = layer(x,capture)
            if capture: maps.append(weights)
        return self.norm(x), torch.stack(maps,dim=1) if capture else None


def numpy_head_probabilities(logits):
    """Softmax over donors; all forbidden returns zeros, not NaNs."""
    maximum = np.max(logits,axis=1,keepdims=True)
    maximum = np.where(np.isfinite(maximum),maximum,0)
    values = np.exp(logits.astype('float64')-maximum)
    total = values.sum(1,keepdims=True)
    return np.divide(values,total,out=np.zeros_like(values),where=total>0)


def numpy_mixture_scores(logits,gate):
    probability = np.sum(numpy_head_probabilities(logits)*gate[:,None,:],axis=-1)
    with np.errstate(divide='ignore'):
        return np.log(probability)


def torch_mixture_scores(logits,gate):
    return torch.logsumexp(torch.log_softmax(logits,dim=1)+gate.clamp_min(1e-30).log()[:,None,:],dim=-1)
