from __future__ import annotations
import numpy as np

def init_chains_in_box(n_chains,chain_length,box_size,bond_length,rng):
    pos=np.zeros((n_chains*chain_length,3),dtype=np.float64)
    for c in range(n_chains):
        start=rng.uniform(0.0,box_size,size=3)
        pos[c*chain_length]=start
        for j in range(1,chain_length):
            v=rng.normal(size=3)
            n=float(np.linalg.norm(v))
            if n==0.0:
                v=np.array([1.0,0.0,0.0]); n=1.0
            pos[c*chain_length+j]=pos[c*chain_length+j-1]+bond_length*(v/n)
    return pos

def apply_pbc(pos,box_size):
    return pos-box_size*np.floor(pos/box_size)

def bonded_pairs(n_chains,chain_length):
    out=[]
    for c in range(n_chains):
        base=c*chain_length
        for j in range(chain_length-1):
            out.append((base+j,base+j+1))
    return out

def chain_index(bead_idx,chain_length):
    return bead_idx//chain_length
