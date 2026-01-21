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
) -> Dict[int, np.ndarray]:
    """Run one independent trajectory and return snapshots at requested sweep times.

    Times are measured in MC sweeps after the quench (t=0).
    """
    sqrt_epsilon = math.sqrt(max(0.0, epsilon))

    # Init
    pos = init_random_walk(params.N, bond_length=bond_length, rng=rng)

    # Pre-equilibrate at high T
    for _ in range(pre_sweeps):
        metropolis_sweep(pos, beta0, step_size, params, eta, sqrt_epsilon, rng)

    # Measurement phase at beta
    want = set(int(t) for t in snapshot_times)
    snapshots: Dict[int, np.ndarray] = {}

    if 0 in want:
        snapshots[0] = pos.copy()

    for sweep in range(1, meas_sweeps + 1):
        metropolis_sweep(pos, beta, step_size, params, eta, sqrt_epsilon, rng)
        if sweep in want:
            snapshots[sweep] = pos.copy()

    missing = want.difference(snapshots.keys())
    if missing:
        raise RuntimeError(f"Missing snapshots for times: {sorted(missing)[:10]} ...")

    return snapshots


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


def compute_Q_D4_for_trajectory(
    snapshots: Dict[int, np.ndarray],
    tw_list: Sequence[int],
    lag_list: Sequence[int],
    pair_i: np.ndarray,
    pair_j: np.ndarray,
    rc2: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return Q and D4 arrays of shape (len(tw_list), len(lag_list))."""
    # Cache distance matrices and contact vectors per snapshot time
    times = sorted(snapshots.keys())
    dist_cache: Dict[int, np.ndarray] = {}
    cont_cache: Dict[int, np.ndarray] = {}

    for t in times:
        ds = dist_sq_matrix(snapshots[t])
        dist_cache[t] = ds
        cont_cache[t] = contact_vector_from_dist_sq(ds, pair_i, pair_j, rc2)

    Q = np.empty((len(tw_list), len(lag_list)), dtype=np.float64)
    D4 = np.empty((len(tw_list), len(lag_list)), dtype=np.float64)

    for a, tw in enumerate(tw_list):
        ds0 = dist_cache[int(tw)]
        c0 = cont_cache[int(tw)]
        for b, lag in enumerate(lag_list):
            t2 = int(tw) + int(lag)
            ds1 = dist_cache[t2]
            c1 = cont_cache[t2]
            D4[a, b] = float(np.mean((ds1 - ds0) ** 2))
            Q[a, b] = float(np.mean(c0 & c1))

    return Q, D4


def log_lags(t_max: int, n_lags: int) -> np.ndarray:
    """Integer log-spaced lags from 1..t_max (unique, sorted)."""
    raw = np.logspace(0.0, math.log10(float(t_max)), num=int(n_lags))
    lags = np.unique(np.clip(np.rint(raw).astype(int), 1, int(t_max)))
    return lags


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

    # CSV header (append mode; write header only if file empty)
    need_header = not os.path.exists(out_csv) or os.path.getsize(out_csv) == 0

    with open(out_csv, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if need_header:
            writer.writerow(
                [
                    "ensemble",
                    "epsilon",
                    "N",
                    "kappa",
                    "pi",
                    "disorder_idx",
                    "tw",
                    "lag",
                    "Q_mean",
                    "chi4",
                    "D4_mean",
                    "D4_std",
                    "n_traj",
                    "n_pairs",
                ]
            )

        for d in range(run.n_disorder):
            # Disorder realization
            rng_dis = np.random.default_rng(base_seed + 100000 * d + 123)
            sigma_seq: Optional[np.ndarray] = None

            if ensemble == "iid":
                eta = generate_eta_iid(N, rng_dis)
            elif ensemble == "correlated":
                eta, sigma_seq = generate_eta_correlated(N, pi=run.pi, kappa=run.kappa, rng=rng_dis)
            else:
                raise ValueError(f"Unknown ensemble: {ensemble}")

            # Trajectory ensemble at fixed disorder
            Q_list: List[np.ndarray] = []
            D4_list: List[np.ndarray] = []

            for tr in range(run.n_traj):
                rng_tr = np.random.default_rng(base_seed + 100000 * d + 1000 * tr + 999)

                snaps = run_aging_trajectory(
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
                )

                Q, D4 = compute_Q_D4_for_trajectory(
                    snaps,
                    tw_list=run.tw_list,
                    lag_list=lag_list,
                    pair_i=pair_i,
                    pair_j=pair_j,
                    rc2=rc2,
                )

                Q_list.append(Q)
                D4_list.append(D4)

            # Aggregate over trajectories -> chi4
            Q_arr = np.stack(Q_list, axis=0)  # (n_traj, n_tw, n_lag)
            D_arr = np.stack(D4_list, axis=0)

            Q_mean = Q_arr.mean(axis=0)
            Q2_mean = (Q_arr ** 2).mean(axis=0)
            chi4 = n_pairs * (Q2_mean - Q_mean ** 2)

            D4_mean = D_arr.mean(axis=0)
            D4_std = D_arr.std(axis=0, ddof=1) if run.n_traj > 1 else np.zeros_like(D4_mean)

            # Write rows
            for a, tw in enumerate(run.tw_list):
                for b, lag in enumerate(lag_list):
                    writer.writerow(
                        [
                            ensemble,
                            float(run.epsilon),
                            int(N),
                            float(run.kappa) if ensemble == "correlated" else 0.0,
                            float(run.pi) if ensemble == "correlated" else 0.0,
                            int(d),
                            int(tw),
                            int(lag),
                            float(Q_mean[a, b]),
                            float(chi4[a, b]),
                            float(D4_mean[a, b]),
                            float(D4_std[a, b]),
                            int(run.n_traj),
                            int(n_pairs),
                        ]
                    )


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
