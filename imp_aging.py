#!/usr/bin/env python3
"""IMP heteropolymer: Metropolis MC + quench/aging + D4/Q/chi4.

This is a research-code starting point for the following locked spec:

Model (3D off-lattice IMP-style heteropolymer):
  H = sum_{i<j} [ h * d_{ij}^2 * 1(j=i+1)
                  + R / d_{ij}^{12}
                  - A / d_{ij}^{6}
                  + sqrt(epsilon) * eta_{ij} / d_{ij}^{6} ]
where eta_{ij} is quenched disorder (mean 0, var 1, symmetric).

Aging protocol:
  - pre-equilibrate at beta0
  - quench to beta=1 at t=0
  - measure two-time observables at waiting times t_w and lags t

Observables:
  D4(t_w,t) = (1/N^2) * sum_{i,j} ( d_{ij}^2(t_w+t) - d_{ij}^2(t_w) )^2
  Q(t_w,t)  = overlap of contact maps (excluding nearest neighbors)
  chi4      = N_pairs * ( <Q^2> - <Q>^2 ) across independent trajectories

Disorder ensembles:
  (A) iid: eta_{ij} ~ N(0,1)
  (B) correlated: eta_{ij} = kappa*sigma_i*sigma_j + sqrt(1-kappa^2)*xi_{ij}
      where sigma is a Markov chain on {+1,-1} with persistence pi.

Defaults target quick test runs. For publishable statistics, increase
n_disorder and n_traj.

Notes:
- We re-center the center of mass after each sweep (energy is translation
  invariant).
- This code prioritizes clarity; N=30 is small enough that O(N^2) pieces are fine.

References (model context):
- Pliszka & Marinari, "On Heteropolymer Shape Dynamics" (1992).
- Irback & Schwarze, "Sequence Dependence of Self-Interacting Random Chains" (1995).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


# -----------------------------
# Model / parameters
# -----------------------------


@dataclass(frozen=True)
class ModelParams:
    N: int = 30
    A: float = 3.8
    R: float = 2.0
    h: float = 1.0


@dataclass(frozen=True)
class RunParams:
    # Thermodynamic parameters
    beta0: float = 0.05
    beta: float = 1.0

    # Disorder + dynamics
    epsilon: float = 6.0
    step_size: float = 0.25

    # Protocol
    pre_sweeps: int = 2000

    # Sampling
    tw_list: Tuple[int, ...] = (0, 100, 300, 1000, 3000, 10000)
    t_max: int = 10000
    n_lags: int = 28

    # Statistics
    n_disorder: int = 2
    n_traj: int = 4

    # Correlated disorder knobs
    kappa: float = 0.7
    pi: float = 0.9


def lj_rmin(params: ModelParams) -> float:
    """Minimum of deterministic LJ part R/r^12 - A/r^6."""
    return (2.0 * params.R / params.A) ** (1.0 / 6.0)


def contact_cutoff(params: ModelParams, factor: float = 1.25) -> float:
    return factor * lj_rmin(params)


# -----------------------------
# Quenched disorder generation
# -----------------------------


def generate_eta_iid(N: int, rng: np.random.Generator) -> np.ndarray:
    """Symmetric iid Gaussian eta_{ij} with mean 0, var 1."""
    eta = np.zeros((N, N), dtype=np.float64)
    iu = np.triu_indices(N, 1)
    vals = rng.normal(loc=0.0, scale=1.0, size=iu[0].shape[0])
    eta[iu] = vals
    eta[(iu[1], iu[0])] = vals
    return eta


def generate_sigma_markov(N: int, pi: float, rng: np.random.Generator) -> np.ndarray:
    """Binary Markov chain sigma_i in {+1,-1} with persistence pi."""
    if not (0.0 <= pi <= 1.0):
        raise ValueError(f"pi must be in [0,1], got {pi}")
    sigma = np.empty(N, dtype=np.int8)
    sigma[0] = 1 if rng.random() < 0.5 else -1
    for i in range(N - 1):
        if rng.random() < pi:
            sigma[i + 1] = sigma[i]
        else:
            sigma[i + 1] = -sigma[i]
    return sigma


def generate_eta_correlated(
    N: int,
    pi: float,
    kappa: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Correlated disorder eta_{ij} and the underlying sigma sequence."""
    if not (0.0 <= kappa <= 1.0):
        raise ValueError(f"kappa must be in [0,1], got {kappa}")

    sigma = generate_sigma_markov(N, pi=pi, rng=rng)
    eta = np.zeros((N, N), dtype=np.float64)

    iu = np.triu_indices(N, 1)
    xi = rng.normal(loc=0.0, scale=1.0, size=iu[0].shape[0])
    sprod = sigma[iu[0]] * sigma[iu[1]]
    vals = kappa * sprod + math.sqrt(max(0.0, 1.0 - kappa * kappa)) * xi

    eta[iu] = vals
    eta[(iu[1], iu[0])] = vals
    return eta, sigma


# -----------------------------
# Geometry helpers
# -----------------------------


def make_contact_pair_indices(N: int) -> Tuple[np.ndarray, np.ndarray]:
    """Pairs (i<j) excluding nearest neighbors (|i-j|>1)."""
    ii: List[int] = []
    jj: List[int] = []
    for i in range(N):
        for j in range(i + 1, N):
            if (j - i) > 1:
                ii.append(i)
                jj.append(j)
    return np.asarray(ii, dtype=np.int32), np.asarray(jj, dtype=np.int32)


def dist_sq_matrix(pos: np.ndarray) -> np.ndarray:
    """All-pairs squared distances (NxN)."""
    diff = pos[:, None, :] - pos[None, :, :]
    return np.einsum("ijk,ijk->ij", diff, diff)


# -----------------------------
# Energy / MC dynamics
# -----------------------------


def delta_energy_move(
    pos: np.ndarray,
    i: int,
    new_pos_i: np.ndarray,
    params: ModelParams,
    eta: np.ndarray,
    sqrt_epsilon: float,
) -> float:
    """Energy change from moving bead i to new_pos_i (vectorized over j!=i)."""
    N = pos.shape[0]

    # Indices of other beads
    mask = np.ones(N, dtype=bool)
    mask[i] = False
    idx = np.nonzero(mask)[0]

    # Distances to others, old/new
    rj = pos[idx]
    old_vec = pos[i] - rj
    new_vec = new_pos_i - rj

    d2_old = np.einsum("ij,ij->i", old_vec, old_vec)
    d2_new = np.einsum("ij,ij->i", new_vec, new_vec)

    # Avoid division-by-zero explosions
    d2_old = np.maximum(d2_old, 1e-12)
    d2_new = np.maximum(d2_new, 1e-12)

    inv2_old = 1.0 / d2_old
    inv2_new = 1.0 / d2_new

    inv6_old = inv2_old * inv2_old * inv2_old
    inv6_new = inv2_new * inv2_new * inv2_new

    inv12_old = inv6_old * inv6_old
    inv12_new = inv6_new * inv6_new

    eta_vec = eta[i, idx]

    e_old = params.R * inv12_old - params.A * inv6_old + sqrt_epsilon * eta_vec * inv6_old
    e_new = params.R * inv12_new - params.A * inv6_new + sqrt_epsilon * eta_vec * inv6_new

    # Harmonic bonds to neighbors (i-1,i) and (i,i+1) if they exist
    neighbor = (idx == (i - 1)) | (idx == (i + 1))
    if np.any(neighbor):
        e_old[neighbor] += params.h * d2_old[neighbor]
        e_new[neighbor] += params.h * d2_new[neighbor]

    return float(np.sum(e_new - e_old))


def compute_energy(pos,params,eta,sqrt_eps):
    N=pos.shape[0]
    dsq=dist_sq_matrix(pos)
    dsq=np.maximum(dsq,1e-12)
    np.fill_diagonal(dsq,1.0)
    inv2=1.0/dsq
    inv6=inv2*inv2*inv2
    inv12=inv6*inv6
    E_lj=params.R*inv12-params.A*inv6+sqrt_eps*eta*inv6
    np.fill_diagonal(E_lj,0.0)
    E_tot=0.5*np.sum(E_lj)
    for i in range(N-1):
        d2=float(np.sum((pos[i+1]-pos[i])**2))
        E_tot+=params.h*d2
    return E_tot


def compute_Rg(pos):
    cm=pos.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((pos-cm)**2,axis=1))))


def compute_ncontacts(pos,pair_i,pair_j,rc2):
    dsq=dist_sq_matrix(pos)
    return int(np.sum(dsq[pair_i,pair_j]<=rc2))


def metropolis_sweep(
    pos: np.ndarray,
    beta: float,
    step_size: float,
    params: ModelParams,
    eta: np.ndarray,
    sqrt_epsilon: float,
    rng: np.random.Generator,
) -> float:
    """One MC sweep (N attempted single-bead moves). Returns acceptance fraction."""
    N = pos.shape[0]
    accepted = 0

    for _ in range(N):
        i = int(rng.integers(0, N))
        delta = rng.uniform(-step_size, step_size, size=3)
        new_pos = pos[i] + delta

        dE = delta_energy_move(pos, i, new_pos, params, eta, sqrt_epsilon)

        if dE <= 0.0:
            pos[i] = new_pos
            accepted += 1
        else:
            # Metropolis accept
            if rng.random() < math.exp(-beta * dE):
                pos[i] = new_pos
                accepted += 1

    # Fix center-of-mass at origin
    pos -= pos.mean(axis=0)
    return accepted / float(N)


def init_random_walk(N: int, bond_length: float, rng: np.random.Generator) -> np.ndarray:
    """Simple random-walk initialization with fixed step length."""
    pos = np.zeros((N, 3), dtype=np.float64)
    for i in range(1, N):
        v = rng.normal(size=3)
        n = float(np.linalg.norm(v))
        if n == 0.0:
            v = np.array([1.0, 0.0, 0.0])
            n = 1.0
        pos[i] = pos[i - 1] + bond_length * (v / n)
    pos -= pos.mean(axis=0)
    return pos


def run_aging_trajectory(
    params: ModelParams,
    eta: np.ndarray,
    epsilon: float,
    beta0: float,
    beta: float,
    pre_sweeps: int,
    meas_sweeps: int,
    step_size: float,
    snapshot_times: Sequence[int],
    rng: np.random.Generator,
    bond_length: float = 1.0,
    pair_i: Optional[np.ndarray] = None,
    pair_j: Optional[np.ndarray] = None,
    rc2: float = 0.0,
) -> Dict:
    sqrt_eps = math.sqrt(max(0.0, epsilon))
    pos = init_random_walk(params.N, bond_length=bond_length, rng=rng)
    for _ in range(pre_sweeps):
        metropolis_sweep(pos, beta0, step_size, params, eta, sqrt_eps, rng)
    want = set(int(t) for t in snapshot_times)
    snapshots: Dict[int, np.ndarray] = {}
    energy: Dict[int, float] = {}
    rg: Dict[int, float] = {}
    nc: Dict[int, int] = {}
    def _record(t, p):
        snapshots[t] = p.copy()
        energy[t] = compute_energy(p, params, eta, sqrt_eps)
        rg[t] = compute_Rg(p)
        if pair_i is not None and pair_j is not None:
            nc[t] = compute_ncontacts(p, pair_i, pair_j, rc2)
    if 0 in want:
        _record(0, pos)
    for sweep in range(1, meas_sweeps + 1):
        metropolis_sweep(pos, beta, step_size, params, eta, sqrt_eps, rng)
        if sweep in want:
            _record(sweep, pos)
    missing = want.difference(snapshots.keys())
    if missing:
        raise RuntimeError(f"Missing snapshots for times: {sorted(missing)[:10]} ...")
    return {"snapshots":snapshots,"energy":energy,"Rg":rg,"n_contacts":nc}


# -----------------------------
# Observables
# -----------------------------


def contact_vector_from_dist_sq(
    dsq: np.ndarray,
    pair_i: np.ndarray,
    pair_j: np.ndarray,
    rc2: float,
) -> np.ndarray:
    return dsq[pair_i, pair_j] <= rc2


def compute_observables_for_trajectory(
    snapshots: Dict[int, np.ndarray],
    tw_list: Sequence[int],
    lag_list: Sequence[int],
    pair_i: np.ndarray,
    pair_j: np.ndarray,
    rc2: float,
) -> Dict[str, np.ndarray]:
    times = sorted(snapshots.keys())
    dist_cache: Dict[int, np.ndarray] = {}
    cont_cache: Dict[int, np.ndarray] = {}
    for t in times:
        ds = dist_sq_matrix(snapshots[t])
        dist_cache[t] = ds
        cont_cache[t] = contact_vector_from_dist_sq(ds, pair_i, pair_j, rc2)
    ntw=len(tw_list)
    nlag=len(lag_list)
    Q=np.empty((ntw,nlag),dtype=np.float64)
    D4=np.empty((ntw,nlag),dtype=np.float64)
    MSD=np.empty((ntw,nlag),dtype=np.float64)
    a2=np.empty((ntw,nlag),dtype=np.float64)
    for a,tw in enumerate(tw_list):
        ds0=dist_cache[int(tw)]
        c0=cont_cache[int(tw)]
        p0=snapshots[int(tw)]
        for b,lag in enumerate(lag_list):
            t2=int(tw)+int(lag)
            ds1=dist_cache[t2]
            c1=cont_cache[t2]
            p1=snapshots[t2]
            D4[a,b]=float(np.mean((ds1-ds0)**2))
            Q[a,b]=float(np.mean(c0&c1))
            disp=p1-p0
            r2=np.sum(disp**2,axis=1)
            MSD[a,b]=float(np.mean(r2))
            r4=r2**2
            mr2=float(np.mean(r2))
            mr4=float(np.mean(r4))
            if mr2>0:
                a2[a,b]=0.6*mr4/(mr2*mr2)-1.0
            else:
                a2[a,b]=0.0
    return {"Q":Q,"D4":D4,"MSD":MSD,"alpha2":a2}


def log_lags(t_max: int, n_lags: int) -> np.ndarray:
    raw = np.logspace(0.0, math.log10(float(t_max)), num=int(n_lags))
    lags = np.unique(np.clip(np.rint(raw).astype(int), 1, int(t_max)))
    return np.concatenate([[0], lags])


def build_snapshot_times(tw_list: Sequence[int], lag_list: Sequence[int]) -> List[int]:
    times = set(int(tw) for tw in tw_list)
    for tw in tw_list:
        for lag in lag_list:
            times.add(int(tw) + int(lag))
    return sorted(times)


# -----------------------------
# Driver
# -----------------------------


def run_condition(
    *,
    params: ModelParams,
    run: RunParams,
    ensemble: str,
    base_seed: int,
    out_csv: str,
) -> None:
    """Run one (epsilon, ensemble) condition and append rows to out_csv."""

    N = params.N
    rc = contact_cutoff(params)
    rc2 = rc * rc
    pair_i, pair_j = make_contact_pair_indices(N)
    n_pairs = int(pair_i.shape[0])

    lag_list = log_lags(run.t_max, run.n_lags)
    snapshot_times = build_snapshot_times(run.tw_list, lag_list)
    meas_sweeps = int(max(snapshot_times))

    # Write a tiny metadata json next to output for reproducibility
    meta_path = os.path.splitext(out_csv)[0] + f"_{ensemble}_eps{run.epsilon:g}.meta.json"
    meta = {
        "model_params": asdict(params),
        "run_params": asdict(run),
        "ensemble": ensemble,
        "rc": rc,
        "n_pairs": n_pairs,
        "lag_list": lag_list.tolist(),
        "snapshot_times": snapshot_times,
        "base_seed": base_seed,
    }
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)

    need_header = not os.path.exists(out_csv) or os.path.getsize(out_csv) == 0
    with open(out_csv, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if need_header:
            writer.writerow([
                "ensemble","epsilon","N","kappa","pi","disorder_idx",
                "tw","lag","Q_mean","chi4","D4_mean","D4_std",
                "n_traj","n_pairs",
                "MSD_mean","alpha2_mean","energy_mean","Rg_mean","n_contacts_mean",
            ])
        for d in range(run.n_disorder):
            rng_dis = np.random.default_rng(base_seed + 100000 * d + 123)
            sigma_seq: Optional[np.ndarray] = None
            if ensemble == "iid":
                eta = generate_eta_iid(N, rng_dis)
            elif ensemble == "correlated":
                eta, sigma_seq = generate_eta_correlated(N, pi=run.pi, kappa=run.kappa, rng=rng_dis)
            else:
                raise ValueError(f"Unknown ensemble: {ensemble}")
            Q_list: List[np.ndarray] = []
            D4_list: List[np.ndarray] = []
            MSD_list: List[np.ndarray] = []
            a2_list: List[np.ndarray] = []
            E_list: List[Dict[int,float]] = []
            Rg_list: List[Dict[int,float]] = []
            nc_list: List[Dict[int,int]] = []
            for tr in range(run.n_traj):
                rng_tr = np.random.default_rng(base_seed + 100000 * d + 1000 * tr + 999)
                traj = run_aging_trajectory(
                    params=params,
                    eta=eta,
                    epsilon=run.epsilon,
                    beta0=run.beta0,
                    beta=run.beta,
                    pre_sweeps=run.pre_sweeps,
                    meas_sweeps=meas_sweeps,
                    step_size=run.step_size,
                    snapshot_times=snapshot_times,
                    rng=rng_tr,
                    bond_length=lj_rmin(params),
                    pair_i=pair_i,
                    pair_j=pair_j,
                    rc2=rc2,
                )
                obs = compute_observables_for_trajectory(
                    traj["snapshots"],
                    tw_list=run.tw_list,
                    lag_list=lag_list,
                    pair_i=pair_i,
                    pair_j=pair_j,
                    rc2=rc2,
                )
                Q_list.append(obs["Q"])
                D4_list.append(obs["D4"])
                MSD_list.append(obs["MSD"])
                a2_list.append(obs["alpha2"])
                E_list.append(traj["energy"])
                Rg_list.append(traj["Rg"])
                nc_list.append(traj["n_contacts"])
            Q_arr = np.stack(Q_list, axis=0)
            D_arr = np.stack(D4_list, axis=0)
            MSD_arr = np.stack(MSD_list, axis=0)
            a2_arr = np.stack(a2_list, axis=0)
            Q_mean = Q_arr.mean(axis=0)
            Q2_mean = (Q_arr**2).mean(axis=0)
            chi4 = n_pairs*(Q2_mean-Q_mean**2)
            D4_mean = D_arr.mean(axis=0)
            D4_std = D_arr.std(axis=0, ddof=1) if run.n_traj > 1 else np.zeros_like(D4_mean)
            MSD_mean = MSD_arr.mean(axis=0)
            a2_mean = a2_arr.mean(axis=0)
            for a, tw in enumerate(run.tw_list):
                tw_int=int(tw)
                for b, lag in enumerate(lag_list):
                    t2=tw_int+int(lag)
                    E_m=float(np.mean([el[t2] for el in E_list]))
                    Rg_m=float(np.mean([rl[t2] for rl in Rg_list]))
                    nc_m=float(np.mean([nl[t2] for nl in nc_list]))
                    writer.writerow([
                        ensemble,
                        float(run.epsilon),
                        int(N),
                        float(run.kappa) if ensemble=="correlated" else 0.0,
                        float(run.pi) if ensemble=="correlated" else 0.0,
                        int(d),
                        int(tw),
                        int(lag),
                        float(Q_mean[a,b]),
                        float(chi4[a,b]),
                        float(D4_mean[a,b]),
                        float(D4_std[a,b]),
                        int(run.n_traj),
                        int(n_pairs),
                        float(MSD_mean[a,b]),
                        float(a2_mean[a,b]),
                        E_m,
                        Rg_m,
                        nc_m,
                    ])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="IMP heteropolymer aging + chi4")

    p.add_argument("--out", type=str, default="/mnt/data/imp_aging_sim/results.csv", help="Output CSV")
    p.add_argument("--seed", type=int, default=12345, help="Base RNG seed")

    # Model parameters
    p.add_argument("--N", type=int, default=30, help="Number of monomers in the chain")

    p.add_argument("--n_disorder", type=int, default=2)
    p.add_argument("--n_traj", type=int, default=4)

    p.add_argument("--epsilon", type=float, nargs="+", default=[0.0, 3.0, 6.0])
    p.add_argument("--ensemble", type=str, nargs="+", default=["iid", "correlated"], choices=["iid", "correlated"])

    p.add_argument("--beta0", type=float, default=0.05)
    p.add_argument("--beta", type=float, default=1.0)
    p.add_argument("--step", type=float, default=0.25)
    p.add_argument("--pre_sweeps", type=int, default=2000)

    p.add_argument("--tw", type=int, nargs="+", default=[0, 100, 300, 1000, 3000, 10000])
    p.add_argument("--t_max", type=int, default=10000)
    p.add_argument("--n_lags", type=int, default=28)

    p.add_argument("--kappa", type=float, default=0.7)
    p.add_argument("--pi", type=float, default=0.9)

    return p.parse_args()


def main() -> None:
    args = parse_args()

    params = ModelParams(N=int(args.N), A=3.8, R=2.0, h=1.0)

    for ens in args.ensemble:
        for eps in args.epsilon:
            run = RunParams(
                beta0=float(args.beta0),
                beta=float(args.beta),
                epsilon=float(eps),
                step_size=float(args.step),
                pre_sweeps=int(args.pre_sweeps),
                tw_list=tuple(int(x) for x in args.tw),
                t_max=int(args.t_max),
                n_lags=int(args.n_lags),
                n_disorder=int(args.n_disorder),
                n_traj=int(args.n_traj),
                kappa=float(args.kappa),
                pi=float(args.pi),
            )

            print(f"Running ensemble={ens} epsilon={eps} -> {args.out}")
            run_condition(
                params=params,
                run=run,
                ensemble=ens,
                base_seed=int(args.seed),
                out_csv=args.out,
            )

    print("Done.")


if __name__ == "__main__":
    main()
