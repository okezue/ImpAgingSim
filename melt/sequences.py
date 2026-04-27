from __future__ import annotations
import numpy as np

def generate_random(N,f_A,rng):
    return (rng.random(N)<f_A).astype(np.int8)

def generate_block(N,block_length,f_A=0.5):
    s=np.empty(N,dtype=np.int8)
    nA=int(round(block_length*f_A*2))
    nA=max(1,nA)
    nB=max(1,2*block_length-nA)
    period=nA+nB
    for i in range(N):
        s[i]=1 if (i%period)<nA else 0
    return s

def generate_alternating(N):
    s=np.zeros(N,dtype=np.int8)
    s[::2]=1
    return s

def generate_correlated(N,kappa,pi,rng):
    if not(0.0<=pi<=1.0):
        raise ValueError(f"pi must be in [0,1], got {pi}")
    if not(0.0<=kappa<=1.0):
        raise ValueError(f"kappa must be in [0,1], got {kappa}")
    sigma=np.empty(N,dtype=np.int8)
    sigma[0]=1 if rng.random()<0.5 else -1
    for i in range(N-1):
        sigma[i+1]=sigma[i] if rng.random()<pi else -sigma[i]
    if kappa<1.0:
        flip=rng.random(N)>kappa
        rand=np.where(rng.random(N)<0.5,1,-1).astype(np.int8)
        sigma=np.where(flip,rand,sigma).astype(np.int8)
    return ((sigma+1)//2).astype(np.int8)

def autocorrelation(seq,kmax):
    s=2*seq.astype(np.float64)-1.0
    N=len(s)
    out=np.empty(kmax+1)
    for k in range(kmax+1):
        out[k]=float(np.mean(s[:N-k]*s[k:])) if k<N else 0.0
    return out
