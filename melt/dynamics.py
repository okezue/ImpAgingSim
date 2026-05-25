"""Two-time dynamic observables for glass-aging analysis.

Standard observables for supercooled/glass-aging studies, computed from saved per-bead
trajectories. Used to support claims about whether sequence correlation suppresses glassy
slow dynamics. Replaces the static-contrast drift used in the previous aging analysis,
which is not by itself a glass-aging diagnostic.

References:
- Berthier & Biroli, Rev. Mod. Phys. 83, 587 (2011).
- Kob & Andersen, Phys. Rev. E 51, 4626 (1995).
- Lacevic, Starr, Schroder, Glotzer, J. Chem. Phys. 119, 7372 (2003) for chi_4.
"""
from __future__ import annotations
import numpy as np

def unwrap_periodic(pos_traj,box_size):
    """Remove PBC jumps assuming consecutive snapshots are close. Returns unwrapped positions
    with the same shape as input (T, N, 3). Used so MSD/F_s are computed on true displacements."""
    L=float(box_size)
    out=np.asarray(pos_traj,dtype=np.float64).copy()
    for t in range(1,out.shape[0]):
        d=out[t]-out[t-1]
        out[t]-=L*np.round(d/L)+0.0
        out[t]=out[t-1]+(out[t]-out[t-1])
    return out

def msd_from_tw(pos_unwrapped,tw_idx):
    """Mean-squared displacement vs lag time, referenced to waiting-time index tw_idx.
    Returns 1D array of length T-tw_idx with MSD(lag) averaged over beads."""
    p=pos_unwrapped
    T=p.shape[0]
    ref=p[tw_idx]
    out=np.zeros(T-tw_idx)
    for s in range(T-tw_idx):
        d=p[tw_idx+s]-ref
        out[s]=float(np.mean(np.sum(d*d,axis=-1)))
    return out

def fs_self_intermediate(pos_unwrapped,tw_idx,k_mag,n_directions=12,rng=None):
    """Self-intermediate scattering function F_s(k, lag; t_w) averaged over n_directions
    random unit vectors of magnitude k_mag and over all beads."""
    if rng is None:
        rng=np.random.default_rng(0)
    p=pos_unwrapped
    T,N,_=p.shape
    ref=p[tw_idx]
    dirs=rng.normal(size=(n_directions,3))
    dirs/=np.linalg.norm(dirs,axis=-1,keepdims=True)
    kvecs=k_mag*dirs
    fs=np.zeros(T-tw_idx)
    for s in range(T-tw_idx):
        d=p[tw_idx+s]-ref
        phases=d@kvecs.T
        fs[s]=float(np.mean(np.cos(phases)))
    return fs

def alpha2_non_gaussian(pos_unwrapped,tw_idx):
    """Non-Gaussian parameter alpha_2(lag; t_w) = 3<r^4>/(5<r^2>^2) - 1. Peaks around
    the alpha-relaxation timescale; zero for Gaussian dynamics."""
    p=pos_unwrapped
    T=p.shape[0]
    ref=p[tw_idx]
    out=np.zeros(T-tw_idx)
    for s in range(T-tw_idx):
        d=p[tw_idx+s]-ref
        r2=np.sum(d*d,axis=-1)
        m2=float(np.mean(r2))
        m4=float(np.mean(r2*r2))
        out[s]=3.0*m4/(5.0*m2*m2)-1.0 if m2>0 else 0.0
    return out

def overlap_self(pos_unwrapped,tw_idx,a=0.3):
    """Self-overlap Q_self(lag; t_w) = (1/N) sum_i w(|r_i(t_w+lag)-r_i(t_w)|), w(r)=Theta(a-r).
    Number fraction of beads still within distance a of their position at t_w. Standard
    glass-aging observable; decays from 1 toward 0 over the alpha-relaxation time."""
    p=pos_unwrapped
    T,N,_=p.shape
    ref=p[tw_idx]
    out=np.zeros(T-tw_idx)
    for s in range(T-tw_idx):
        d=p[tw_idx+s]-ref
        r=np.sqrt(np.sum(d*d,axis=-1))
        out[s]=float(np.mean(r<a))
    return out

def chi4_from_overlap_realizations(overlap_runs,n_particles):
    """chi_4(lag) = N * Var_runs[Q_self(lag)] with N = number of particles in each realization.

    A previous version of this function used N = n_runs (number of independent
    realizations) instead of the system size, which underestimated chi_4 by a factor
    n_runs / N. For the corrected aging campaign (N=5760, n_runs=3) the underestimate
    was a factor of ~1920.

    Parameters
    ----------
    overlap_runs : array_like
        Shape (n_runs, n_lags). Q_self curves from independent realizations.
    n_particles : int
        Number of particles summed in each Q_self.

    Returns
    -------
    ndarray of shape (n_lags,) or nan-filled if n_runs<2.
    """
    Q=np.asarray(overlap_runs,dtype=float)
    if Q.ndim!=2 or Q.shape[0]<2:
        return np.full(Q.shape[-1] if Q.ndim>=1 else 0,np.nan)
    if n_particles<=0:
        raise ValueError("n_particles must be positive")
    return float(n_particles)*np.var(Q,axis=0,ddof=1)

def tau_alpha_from_overlap(lags,Q,threshold=1.0/np.e):
    """Extract alpha-relaxation time from a self-overlap curve as the lag at which Q crosses
    threshold (default 1/e). Returns NaN if the threshold is not crossed within the window."""
    Q=np.asarray(Q)
    lags=np.asarray(lags)
    below=np.where(Q<threshold)[0]
    if len(below)==0:
        return float("nan")
    i=int(below[0])
    if i==0:
        return float(lags[0])
    x0,x1=float(lags[i-1]),float(lags[i])
    y0,y1=float(Q[i-1]),float(Q[i])
    if y0==y1:
        return x1
    return x0+(threshold-y0)*(x1-x0)/(y1-y0)

def kww_fit(lags,Q,p0=None):
    """Fit Q(lag) = A * exp(-(lag/tau)^beta) to the relaxation tail. Returns (A, tau, beta)
    or NaN tuple on failure. Uses scipy.optimize.curve_fit if available, otherwise log-log
    fallback."""
    try:
        from scipy.optimize import curve_fit
    except Exception:
        return (float("nan"),)*3
    lags=np.asarray(lags,dtype=np.float64)
    Q=np.asarray(Q,dtype=np.float64)
    mask=(Q>1e-3)&(Q<1.0)&np.isfinite(Q)&(lags>0)
    if mask.sum()<4:
        return (float("nan"),)*3
    def model(t,A,tau,beta):
        return A*np.exp(-(t/tau)**beta)
    if p0 is None:
        p0=[float(Q[mask][0]),float(lags[mask][len(lags[mask])//2]),0.7]
    try:
        popt,_=curve_fit(model,lags[mask],Q[mask],p0=p0,maxfev=5000)
        return tuple(float(x) for x in popt)
    except Exception:
        return (float("nan"),)*3

def compute_all_dynamics(pos_traj,box_size,tw_indices,k_star,n_directions=12,a=0.3,rng=None):
    """Wrapper computing MSD, F_s(k_star), alpha_2, Q_self for each waiting-time index."""
    if rng is None:
        rng=np.random.default_rng(0)
    pos=unwrap_periodic(pos_traj,box_size)
    out={}
    for tw in tw_indices:
        out[int(tw)]={
            "msd":msd_from_tw(pos,tw),
            "fs":fs_self_intermediate(pos,tw,k_star,n_directions=n_directions,rng=rng),
            "alpha2":alpha2_non_gaussian(pos,tw),
            "Q_self":overlap_self(pos,tw,a=a),
        }
    return out
