"""Explicit TRACE CUDA backend; keep upstream initialization, recursion and ECM."""
import ctypes as ct
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from scipy.stats import gamma as gamma_sp

_lib=None

def library():
    global _lib
    if _lib is not None:
        return _lib
    home=Path(os.environ.get('TRACE_GPU_HOME','/mnt/d/software/gu/trace-gpu'))
    import tracehmm.trace as upstream
    import trace_utils
    manifest=json.loads((home/'build-manifest.json').read_text())
    # An installed TRACE update must be revalidated, not silently monkeypatched.
    for name,digest in manifest['sha256'].items():
        p=Path(name)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
            raise RuntimeError(f'TRACE GPU build/source mismatch: {p}; rebuild and validate')
    lib=ct.CDLL(str(home/'libtrace_gpu.so'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    lib.gu_trace_error.restype=ct.c_char_p
    lib.gu_trace_check.argtypes=[ct.c_int]
    lib.gu_trace_emissions.argtypes=[ct.c_int,arr,arr,arr]+[ct.c_double]*12
    lib.gu_trace_transitions.argtypes=[ct.c_int,arr,arr,arr]+[ct.c_double]*4+[arr,arr]
    def checked(rc):
        if rc:
            raise RuntimeError('TRACE CUDA: '+lib.gu_trace_error().decode())
    lib.checked=checked
    checked(lib.gu_trace_check(int(os.environ.get('TRACE_GPU_DEVICE','0'))))
    _lib=lib
    return lib

def cache_emissions(self,z=1,alpha=1e-2):
    lib=library()
    x=np.ascontiguousarray(self.ncoal,dtype=np.float64)
    mask=np.ascontiguousarray(self.mask,dtype=np.float64)
    if x.ndim!=1 or x.size==0 or mask.shape!=x.shape:
        raise ValueError('Invalid TRACE emission dimensions')
    a1,b1=self.emi2_a1,self.emi2_b1
    a2,b2=(self.emi2_a2,self.emi2_b2) if z==1 else (a1,b1)
    cp=-1.; left=right=0.
    if z==1:
        f=lambda t: gamma_sp.pdf(t,a=a1,scale=1/b1)
        g=lambda t: gamma_sp.pdf(t,a=a2,scale=1/b2)
        _,_,_,crossing=self.find_crossing_point(f=f,g=g,mu1=a2/b2)
        if crossing is not None:
            cp=float(crossing)
            interval=np.max(x)-cp
            gl=gamma_sp.cdf(cp,a=a2,scale=1/b2)
            fl=gamma_sp.cdf(cp,a=a1,scale=1/b1)
            left=float(np.log((fl-alpha)/gl))
            right=float(np.log(alpha/interval))
    out=np.empty_like(x)
    lib.checked(lib.gu_trace_emissions(x.size,x,mask,out,
        a1,b1,a1*math.log(b1),math.lgamma(a1),a2,b2,a2*math.log(b2),math.lgamma(a2),
        cp,left,right,float(np.log(self.scale_factor))))
    if not np.isfinite(out).all():
        raise FloatingPointError('Non-finite TRACE GPU emissions')
    return out

def update_oneind_cython(alphas,betas,emissions,p,q):
    lib=library()
    a=np.ascontiguousarray(alphas,dtype=np.float64)
    b=np.ascontiguousarray(betas,dtype=np.float64)
    if a.ndim!=2 or a.shape[0]!=2 or a.shape!=b.shape or a.shape[1]<1:
        raise ValueError('Invalid TRACE transition dimensions')
    n=a.shape[1]
    # Match the original function's first-n-column indexing, including masked runs.
    e=np.ascontiguousarray(emissions[:,:n],dtype=np.float64)
    if e.shape!=a.shape or not (0<p<1 and 0<q<1):
        raise ValueError('Invalid TRACE transition probabilities/emissions')
    o1=np.empty(n-1); o2=np.empty(n-1)
    lib.checked(lib.gu_trace_transitions(n,a,b,e,math.log(p),math.log(q),math.log(1-p),math.log(1-q),o1,o2))
    return o1,o2

def install():
    library()
    import tracehmm.trace as upstream
    upstream.TRACE.cache_emissions=cache_emissions
    upstream.update_oneind_cython=update_oneind_cython
