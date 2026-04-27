from __future__ import annotations
import os,glob
import numpy as np

def load_run_sf(run_dir):
    p=os.path.join(run_dir,"structure_factor.npz")
    if not os.path.exists(p):
        return None
    d=np.load(p)
    return {"steps":d["steps"].astype(np.int64),"k":d["k"],
            "S_AA":d["S_AA"],"S_BB":d["S_BB"],"S_AB":d["S_AB"]}

def collect_phi_at_kstar(scan_dir,k_star_idx=None):
    runs=[]
    seqs=[]
    seeds=[]
    common_steps=None
    common_k=None
    for run_dir in sorted(glob.glob(os.path.join(scan_dir,"*"))):
        if not os.path.isdir(run_dir):continue
        d=load_run_sf(run_dir)
        if d is None:continue
        if common_steps is None:
            common_steps=d["steps"];common_k=d["k"]
        elif len(d["steps"])!=len(common_steps):
            continue
        runs.append(d["S_AA"])
        bn=os.path.basename(run_dir)
        seq=bn.rsplit("_",1)[0];sd=int(bn.rsplit("_s",1)[-1])
        seqs.append(seq);seeds.append(sd)
    if not runs:
        return None
    S=np.stack(runs,axis=0)
    if k_star_idx is None:
        avg_last=S[:,-3:,:].mean(axis=(0,1))
        k_min_floor=2
        k_star_idx=int(np.argmax(avg_last[k_min_floor:]))+k_min_floor
    return {"steps":common_steps,"k":common_k,"k_star_idx":k_star_idx,
            "S_AA":S,"sequences":seqs,"seeds":seeds}

def two_time_overlap(S_traj,k_star_idx,t_w_indices):
    n_traj,n_t,n_k=S_traj.shape
    A=S_traj[:,:,k_star_idx]
    A=A-A.mean(axis=0,keepdims=True)
    A=A/(np.std(A,axis=0,keepdims=True)+1e-12)
    Q=np.zeros((len(t_w_indices),n_t),dtype=np.float64)
    for a,iw in enumerate(t_w_indices):
        if iw>=n_t:continue
        ref=A[:,iw]
        for b in range(n_t-iw):
            Q[a,b]=float(np.mean(ref*A[:,iw+b]))
    return Q

def chi4_from_traj(S_traj,k_star_idx,t_w_indices):
    n_traj,n_t,n_k=S_traj.shape
    A=S_traj[:,:,k_star_idx]
    chi4=np.zeros((len(t_w_indices),n_t),dtype=np.float64)
    for a,iw in enumerate(t_w_indices):
        if iw>=n_t:continue
        ref=A[:,iw]
        for b in range(n_t-iw):
            prod=ref*A[:,iw+b]
            chi4[a,b]=float(np.var(prod,ddof=1)) if n_traj>1 else 0.0
    return chi4

def compute_twotime_for_scan(scan_dir,t_w_fractions=(0.0,0.25,0.5,0.75)):
    seqs_uniq=set()
    grouped={}
    for run_dir in sorted(glob.glob(os.path.join(scan_dir,"*"))):
        if not os.path.isdir(run_dir):continue
        d=load_run_sf(run_dir)
        if d is None:continue
        bn=os.path.basename(run_dir)
        seq=bn.rsplit("_s",1)[0]
        grouped.setdefault(seq,[]).append((d["steps"],d["S_AA"],d["k"]))
    out={}
    for seq,trajs in grouped.items():
        steps=trajs[0][0];k=trajs[0][2]
        n_t=len(steps);n_k=len(k)
        S=np.zeros((len(trajs),n_t,n_k),dtype=np.float64)
        for i,(_,SAA,_) in enumerate(trajs):
            if SAA.shape==(n_t,n_k):S[i]=SAA
        avg_last=S[:,-3:,:].mean(axis=(0,1))
        kstar_idx=int(np.argmax(avg_last[2:]))+2
        tw_idx=[int(round(f*(n_t-1))) for f in t_w_fractions]
        Q=two_time_overlap(S,kstar_idx,tw_idx)
        chi4=chi4_from_traj(S,kstar_idx,tw_idx)
        out[seq]={"steps":steps,"k":k,"kstar_idx":kstar_idx,
                  "k_star":float(k[kstar_idx]),"tw_idx":tw_idx,
                  "Q":Q,"chi4":chi4,"n_traj":len(trajs)}
    return out
