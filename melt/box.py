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

def relax_overlaps(pos,box_size,chain_length,r_min=0.8,iterations=200,step=0.5):
    """Push apart non-bonded beads closer than r_min, under the minimum image.

    init_chains_in_box lays down independent random walks with no overlap
    rejection, so beads can land essentially coincident.  The WCA core goes as
    r^-12, so a 0.007 sigma contact carries ~1e25 energy and the local energy
    minimizer cannot always recover: such runs diverge within the first
    thousand steps.  This performs a deterministic, RNG-free geometric
    push-off before minimization.  It is a repair of the initial condition
    only; the subsequent high-temperature equilibration sets the ensemble.

    Returns (positions, info) with the minimum non-bonded separation before
    and after, and whether the target was reached.
    """
    from scipy.spatial import cKDTree
    L=float(box_size)
    p=np.mod(np.asarray(pos,dtype=np.float64).copy(),L)
    def min_sep(w):
        t=cKDTree(w,boxsize=L)
        pr=t.query_pairs(r_min,output_type='ndarray')
        if len(pr)==0:
            return np.inf,pr
        keep=~((pr[:,1]-pr[:,0]==1)&(pr[:,0]//chain_length==pr[:,1]//chain_length))
        pr=pr[keep]
        if len(pr)==0:
            return np.inf,pr
        d=(p[pr[:,0]]-p[pr[:,1]]+L/2)%L-L/2
        return float(np.linalg.norm(d,axis=1).min()),pr
    before,_=min_sep(p)
    for _ in range(iterations):
        cur,pr=min_sep(p)
        if not np.isfinite(cur) or cur>=r_min:
            break
        d=(p[pr[:,0]]-p[pr[:,1]]+L/2)%L-L/2
        dist=np.linalg.norm(d,axis=1)
        safe=np.where(dist>1e-9,dist,1e-9)
        unit=d/safe[:,None]
        # coincident beads have no defined direction; push along x
        unit[dist<=1e-9]=np.array([1.0,0.0,0.0])
        push=(step*(r_min-safe))[:,None]*unit
        np.add.at(p,pr[:,0], push/2)
        np.add.at(p,pr[:,1],-push/2)
        p=np.mod(p,L)
    after,_=min_sep(p)
    return p,{"min_separation_before":before,"min_separation_after":after,
              "target_r_min":float(r_min),
              "reached_target":bool(after>=0.99*r_min)}
