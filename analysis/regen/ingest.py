from __future__ import annotations
import json,re
from pathlib import Path
from typing import Iterable
import numpy as np
import pandas as pd

def load_run(rd:Path):
    mp=rd/"meta.json";sf=rd/"structure_factor.npz";sn=rd/"snapshots.csv"
    if not(mp.exists() and sf.exists()):return None
    m=json.loads(mp.read_text())
    s=np.load(sf,allow_pickle=False)
    snap=pd.read_csv(sn) if sn.exists() else None
    return {"meta":m,"sf":{k:s[k] for k in s.files},"snap":snap,"dir":rd}

def walk_chunk(chunk_dir:Path):
    out=[]
    if not chunk_dir.exists():return out
    for rd in sorted(chunk_dir.glob("**/")):
        if rd==chunk_dir:continue
        r=load_run(rd)
        if r is not None:out.append(r)
    return out

def chunk_dataframe(runs:Iterable[dict]):
    rows=[]
    for r in runs:
        m=r["meta"]
        mp=m.get("melt_params",{})
        rp=m.get("run_params",{})
        sf=r["sf"]
        kvec=sf.get("k");sAA=sf.get("S_AA");sAB=sf.get("S_AB");sBB=sf.get("S_BB");steps=sf.get("steps")
        if kvec is None or sAA is None:continue
        n_t=len(steps);last_k=min(5,n_t)
        if n_t==0:continue
        SAA=sAA[-last_k:].mean(0);SAB=sAB[-last_k:].mean(0) if sAB is not None else np.zeros_like(SAA)
        L=float(mp.get("box_size",0.0))
        k_min=1.5*(2*np.pi/L) if L>0 else 0.0
        mask=(kvec>=k_min)&np.isfinite(SAA)&np.isfinite(SAB)
        if not mask.any():continue
        C=SAA-SAB
        i=int(np.argmax(C[mask]))
        k_star=float(kvec[mask][i]);C_star=float(C[mask][i]);S_star=float(SAA[mask][i])
        snap=r["snap"]
        pe=float(snap["total_energy"].tail(5).mean()) if snap is not None and "total_energy" in snap else np.nan
        rg=float(snap["mean_Rg"].tail(5).mean()) if snap is not None and "mean_Rg" in snap else np.nan
        Tinst=float(snap["temperature_inst"].tail(5).mean()) if snap is not None and "temperature_inst" in snap else np.nan
        rows.append({
            "run":r["dir"].name,
            "sequence":m.get("sequence_type",""),
            "n_chains":int(mp.get("n_chains",0)),
            "N":int(mp.get("chain_length",0)),
            "L":L,
            "rho":int(mp.get("n_chains",0))*int(mp.get("chain_length",0))/(L**3) if L>0 else np.nan,
            "kappa":_first_float(snap,"kappa") if snap is not None else np.nan,
            "pi":_first_float(snap,"pi") if snap is not None else np.nan,
            "f_A":_first_float(snap,"f_A") if snap is not None else np.nan,
            "eps_AA":float(mp.get("lj_eps_AA",1.0)),
            "eps_BB":float(mp.get("lj_eps_BB",1.0)),
            "eps_AB":float(mp.get("lj_eps_AB",0.1)),
            "T_quench_K":float(m.get("T_quench_kelvin",mp.get("temperature",np.nan))),
            "T_quench_star":float(m.get("T_quench_star",np.nan)),
            "C":C_star,"k_star":k_star,"xi":2*np.pi/k_star if k_star>0 else np.nan,
            "S_AA_peak":S_star,"E":pe,"Rg":rg,"T_inst":Tinst,
            "seed":int(rp.get("seed",-1)),
        })
    return pd.DataFrame(rows)

def _first_float(snap,col):
    if snap is None or col not in snap:return np.nan
    try:return float(snap[col].iloc[0])
    except Exception:return np.nan

def aggregate(df:pd.DataFrame,group_cols:list[str],value_cols:list[str]):
    g=df.groupby(group_cols,dropna=False)
    out=g[value_cols].agg(["mean","sem","count"]).reset_index()
    out.columns=["_".join([c for c in t if c]).rstrip("_") for t in out.columns]
    return out

def bootstrap_ci(x,n=2000,p=0.95,rng=None):
    if rng is None:rng=np.random.default_rng(0)
    x=np.asarray(x,dtype=np.float64)
    x=x[np.isfinite(x)]
    if len(x)<2:return (np.nan,np.nan,np.nan)
    m=np.array([rng.choice(x,size=len(x),replace=True).mean() for _ in range(n)])
    lo=float(np.quantile(m,(1-p)/2));hi=float(np.quantile(m,1-(1-p)/2))
    return (float(x.mean()),lo,hi)
