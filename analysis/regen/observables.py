from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

def rpa_predictor_finite(kappa,pi,N,k=0.0,b=1.0):
    """Corrected finite-chain bracketed RPA predictor for the bare composition factor:
    1 + 2 kappa^2 * sum_{l=1}^{N-1} (1 - l/N) * (lambda * exp(-k^2 b^2/6))^l, with lambda=2*pi-1.

    At k=0 this is the long-wavelength limit. The lag-zero (1) term is the random-sequence
    baseline 4 f_A(1-f_A) divided out; the kappa^2-weighted geometric tail is the sequence
    memory excess. This is the predictor the reviewer recommended (R^2 jumps from ~0.60
    using the asymptotic kappa^2/[2(1-pi)] alone to ~0.83 using this bracketed quantity)."""
    if not(0.0<=kappa<=1.0):raise ValueError("kappa in [0,1]")
    if not(0.0<=pi<=1.0):raise ValueError("pi in [0,1]")
    lam=2.0*pi-1.0
    if abs(lam)<1e-15:return 1.0
    decay=lam*np.exp(-k*k*b*b/6.0)
    ell=np.arange(1,N)
    w=1.0-ell/N
    tail=2.0*kappa*kappa*np.sum(w*decay**ell)
    return 1.0+tail

def rpa_predictor_asymptotic(kappa,pi):
    """Old asymptotic predictor: kappa^2 / [2*(1-pi)]. Symmetric f_A=1/2 limit only."""
    if pi>=1.0-1e-12:return np.inf
    return kappa*kappa/(2.0*(1.0-pi))

def fit_loglog(x,y,n_boot=2000,rng=None):
    """log10(y) = slope * log10(x) + b. Returns (slope, b, R^2, slope_lo, slope_hi)
    with 95% bootstrap CI on slope. NaNs and non-positive entries are dropped."""
    if rng is None:rng=np.random.default_rng(0)
    x=np.asarray(x,dtype=np.float64);y=np.asarray(y,dtype=np.float64)
    m=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
    if m.sum()<3:return (np.nan,)*5
    lx=np.log10(x[m]);ly=np.log10(y[m])
    s,b=np.polyfit(lx,ly,1)
    yhat=s*lx+b
    ss_res=float(np.sum((ly-yhat)**2));ss_tot=float(np.sum((ly-ly.mean())**2))
    r2=1.0-ss_res/ss_tot if ss_tot>0 else np.nan
    n=len(lx);slopes=np.empty(n_boot)
    for i in range(n_boot):
        idx=rng.integers(0,n,n)
        slopes[i],_=np.polyfit(lx[idx],ly[idx],1)
    lo,hi=np.quantile(slopes,[0.025,0.975])
    return (float(s),float(b),float(r2),float(lo),float(hi))

def collapse_predictor_comparison(df:pd.DataFrame,N=40,b=1.0):
    """For each (pi, kappa) condition in df, attach corrected and asymptotic predictors
    and the excess contrast C - C(kappa=0, same pi, same other params).
    Returns dataframe ready for the Fig 2B side-by-side comparison plot."""
    d=df.copy()
    d=d.dropna(subset=["pi","kappa"])
    base=d[d["kappa"]==0.0].groupby(list(set(d.columns)-{"kappa","C","run","seed","C_lo","C_hi","C_sem"}-{c for c in d.columns if c.startswith("_")}))
    rows=[]
    for keys,grp_kappa_pi in d.groupby([c for c in ["pi","f_A","eps_AB","n_chains","N","L"] if c in d.columns]):
        if isinstance(keys,tuple):key_pi=dict(zip([c for c in ["pi","f_A","eps_AB","n_chains","N","L"] if c in d.columns],keys))
        else:key_pi={"pi":keys}
        c0_rows=grp_kappa_pi[grp_kappa_pi["kappa"]==0.0]
        if c0_rows.empty:continue
        C0=float(c0_rows["C"].mean())
        for _,row in grp_kappa_pi.iterrows():
            k=row["kappa"];p=row["pi"];C=row["C"]
            if not(np.isfinite(k) and np.isfinite(p) and np.isfinite(C)) or k==0.0:continue
            pred_corr=rpa_predictor_finite(k,p,N=N,k=0.0,b=b)
            pred_asy=rpa_predictor_asymptotic(k,p)
            rows.append({**key_pi,"kappa":k,"C":C,"C0":C0,"excess":C-C0,
                         "pred_corrected_minus_1":pred_corr-1.0,"pred_asymptotic":pred_asy})
    return pd.DataFrame(rows)

def aging_slope_bootstrap(tw,C,n_boot=2000,rng=None):
    """Linear regression slope of C vs log10(tw) with bootstrap CI on slope."""
    if rng is None:rng=np.random.default_rng(0)
    tw=np.asarray(tw,dtype=np.float64);C=np.asarray(C,dtype=np.float64)
    m=np.isfinite(tw)&np.isfinite(C)&(tw>0)
    if m.sum()<2:return (np.nan,np.nan,np.nan)
    x=np.log10(tw[m]);y=C[m]
    s,_=np.polyfit(x,y,1)
    n=len(x);ss=np.empty(n_boot)
    for i in range(n_boot):
        idx=rng.integers(0,n,n)
        ss[i],_=np.polyfit(x[idx],y[idx],1)
    lo,hi=np.quantile(ss,[0.025,0.975])
    return (float(s),float(lo),float(hi))

def dynamics_summary(traj_npz_path:Path,box_size,k_star,tw_indices,a=0.3,n_dirs=12,rng=None):
    """Load saved trajectory and compute F_s, MSD, alpha_2, Q_self per waiting-time index.
    Returns dict keyed by tw with each value a dict {msd, fs, alpha2, Q_self, lags}.
    lags are in snapshot-index units; caller converts to BD steps via meta.snapshot_interval."""
    from melt.dynamics import compute_all_dynamics
    d=np.load(traj_npz_path,allow_pickle=False)
    pos=d["positions"]
    if rng is None:rng=np.random.default_rng(0)
    out=compute_all_dynamics(pos,float(box_size),tw_indices,float(k_star),n_directions=n_dirs,a=a,rng=rng)
    for tw,obs in out.items():
        n=len(obs["msd"])
        obs["lags"]=np.arange(n)
    return out

def tau_alpha_per_tw(dyn_out,steps_per_lag,threshold=1.0/np.e):
    """For each tw -> threshold-crossing of F_s gives tau_alpha in BD step units."""
    from melt.dynamics import tau_alpha_from_overlap
    out={}
    for tw,obs in dyn_out.items():
        lags=obs["lags"]*float(steps_per_lag)
        fs=obs["fs"]
        out[tw]=tau_alpha_from_overlap(lags,fs,threshold=threshold)
    return out

def chi4_per_tw(per_seed_Q_arrays):
    """For each tw, stack the self-overlap arrays from independent seeds and compute
    chi_4(t) = n_seeds * Var(Q(t)). Input is dict tw -> 2D array (n_seeds, n_lags)."""
    from melt.dynamics import chi4_from_overlap_realizations
    return {tw:chi4_from_overlap_realizations(arr) for tw,arr in per_seed_Q_arrays.items()}
