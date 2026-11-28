"""Matched linear fit on all regional and global features using bounded tall-matrix batches."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.linalg.blas import dsyrk

RIDGE_PENALTIES = (1., 10., 100., 1000.)


def regional_ridge(array, ix, side, y, build, tune, mu, sd, batch=512):
    # Normal equations accumulated in double precision. No participants-by-participants matrix.
    width = len(mu) + side.shape[1] + 1
    gram = np.zeros((width,width),dtype=float,order="F")
    target = np.zeros(width)

    def batches(rows):
        for start in range(0,len(rows),batch):
            r = rows[start:start+batch]
            x = np.empty((len(r),width),dtype=float)
            x[:,:len(mu)] = ((array[ix[r]]-mu)/sd)[:,:,0]
            x[:,len(mu):-1] = side[r]
            x[:,-1] = 1.
            yield r,x

    for r,x in batches(np.flatnonzero(build)):
        gram = dsyrk(1.,x.T,beta=1.,c=gram,lower=1,overwrite_c=1)
        target += x.T @ y[r]
    coefficients=[]
    for penalty in RIDGE_PENALTIES:
        system=gram.copy(order="F")
        system.ravel(order="F")[::width+1][:-1] += penalty
        coefficients.append(cho_solve(cho_factor(system,lower=True,overwrite_a=True,check_finite=False),
                                      target,check_finite=False))
    candidates=np.column_stack(coefficients)
    sse=np.zeros(len(coefficients))
    for r,x in batches(np.flatnonzero(tune)):
        sse+=np.square(x@candidates-y[r,None]).sum(0)
    best=int(np.argmin(sse))
    beta=coefficients[best]
    raw=np.empty(len(y))
    for r,x in batches(np.arange(len(y))):
        raw[r]=x@beta
    return raw,dict(coefficients=beta,alpha=RIDGE_PENALTIES[best],
                    tune_mse=(sse/tune.sum()).tolist(),penalties=list(RIDGE_PENALTIES))
