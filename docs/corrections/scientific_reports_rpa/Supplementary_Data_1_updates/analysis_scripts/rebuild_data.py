#!/usr/bin/env python3
"""Rebuild the numerical source tables for the manuscript and Supplementary Information.

This script uses only:

1. the 789 processed corrected-run rows supplied with the manuscript archive;
2. the exact symmetric composition reanalysis table supplied in that archive;
3. three recovered corrected trajectories, used only for explicitly labelled
   representative single-seed spectra; and
4. a small set of numerical sensitivity results retained from the available
   processed inputs. Rows from (4) carry an explicit provenance label because the complete
   raw spectrum archive is not present in this workspace.

No potential-energy, density-scaling, waiting-time, aging, or bulk-transition
claim is computed here.
"""

from __future__ import annotations

import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.stats import binomtest, friedmanchisquare, linregress, rankdata


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUT = ROOT / "work" / "heteropolymer_paper_complete" / "source_data"
CORRECTED = INPUT / "source_data_corrected"
TRAJ = (
    ROOT
    / "library_recovered"
    / "source_v4"
    / "nature_physics_heteropolymer_v4_source_data"
    / "selected_trajectories"
)
OUT = HERE / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
RAW_WARNING = OUT / "raw_processed_with_warning"
RAW_WARNING.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "repo"))
from melt.sequences import generate_per_chain  # noqa: E402

G = 56
N_CHAIN = 40
N_CHAINS = 144
N_TOTAL = N_CHAIN * N_CHAINS
L_BOX = 22.0
PI_SELECTED = 0.99
F_A_TARGET = 0.5
K_FUND = 2.0 * np.pi / L_BOX
SCALE_DEFAULT = G**6 / (N_TOTAL * L_BOX**3)


def save(df: pd.DataFrame, name: str) -> None:
    """Write a stable, high-precision CSV."""
    df.to_csv(OUT / name, index=False, float_format="%.15g")


def sem(x: pd.Series | np.ndarray) -> float:
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size < 2:
        return float("nan")
    return float(np.std(a, ddof=1) / np.sqrt(a.size))


def blocked_friedman_permutation(
    values: np.ndarray,
    *,
    n_draws: int = 1_000_000,
    seed: int = 20260803,
    chunk_size: int = 20_000,
) -> tuple[float, float, int]:
    """Monte Carlo Friedman null by independent relabelling within each block.

    ``Generator.permuted(..., axis=2)`` independently permutes every seed row in
    every draw.  Fixed draw count, seed and chunk size make the audit bitwise
    reproducible with the documented NumPy implementation.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not np.all(np.isfinite(values)):
        raise ValueError("values must be a finite blocks-by-conditions matrix")
    n_blocks, n_conditions = values.shape
    ranks = np.apply_along_axis(rankdata, 1, values)
    coefficient = 12.0 / (n_blocks * n_conditions * (n_conditions + 1.0))
    constant = 3.0 * n_blocks * (n_conditions + 1.0)
    observed = coefficient * np.sum(ranks.sum(axis=0) ** 2) - constant

    rng = np.random.default_rng(seed)
    exceedances = 0
    for start in range(0, n_draws, chunk_size):
        n_chunk = min(chunk_size, n_draws - start)
        tiled = np.broadcast_to(ranks, (n_chunk, n_blocks, n_conditions)).copy()
        permuted = rng.permuted(tiled, axis=2)
        statistics = coefficient * np.sum(permuted.sum(axis=1) ** 2, axis=1) - constant
        exceedances += int(np.count_nonzero(statistics >= observed - 1e-12))
    p_value = exceedances / n_draws
    mcse = math.sqrt(p_value * (1.0 - p_value) / n_draws)
    return float(p_value), float(mcse), exceedances


def force_equivalent_kernel_q0() -> tuple[float, float]:
    """Direct spherical quadrature of the unit-depth continuous force kernel."""
    sigma = 1.0
    r_minimum = 2.0 ** (1.0 / 6.0) * sigma
    r_cutoff = 2.5 * sigma

    def u_lj(r: float) -> float:
        sr6 = (sigma / r) ** 6
        return 4.0 * (sr6**2 - sr6)

    cutoff_value = u_lj(r_cutoff)
    inner_value = u_lj(r_minimum) - cutoff_value
    inner_integral, inner_error = quad(
        lambda r: r**2 * inner_value,
        0.0,
        r_minimum,
        epsabs=1e-13,
        epsrel=1e-13,
    )
    tail_integral, tail_error = quad(
        lambda r: r**2 * (u_lj(r) - cutoff_value),
        r_minimum,
        r_cutoff,
        epsabs=1e-13,
        epsrel=1e-13,
    )
    transform_q0 = 4.0 * math.pi * (inner_integral + tail_integral)
    quadrature_error = 4.0 * math.pi * (inner_error + tail_error)
    return float(transform_q0), float(quadrature_error)


def archive_to_number_scale(row: pd.Series) -> float:
    n_tot = int(row["n_chains"]) * int(row["N"])
    return float(G**6 / (n_tot * float(row["L"]) ** 3))


def read_corrected(stem: str) -> pd.DataFrame:
    return pd.read_csv(CORRECTED / f"{stem}.csv")


def sequence_tables() -> None:
    """Reconstruct the exact static-scan sequence inputs from repo code and seeds.

    Production calls seed ``default_rng(seed)`` immediately before
    ``generate_per_chain``.  We reconstruct five static-scan seeds for kappa in
    {0, 0.6, 1}; correlations are averaged over valid within-chain pairs only and
    centered about the known ensemble mean E[psi]=0.
    """
    strip_rows: list[dict] = []
    cov_rows: list[dict] = []
    boundary_rows: list[dict] = []
    for kappa in (0.0, 0.6, 1.0):
        for seed in range(1, 6):
            rng = np.random.default_rng(seed)
            types = generate_per_chain(
                "correlated", N_CHAINS, N_CHAIN, F_A_TARGET, 4, kappa, PI_SELECTED, rng
            ).reshape(N_CHAINS, N_CHAIN)
            psi = 2.0 * types - 1.0
            if seed == 1:
                for chain in range(4):
                    for bead in range(N_CHAIN):
                        strip_rows.append(
                            {
                                "kappa": kappa,
                                "pi": PI_SELECTED,
                                "seed": seed,
                                "chain_index": chain,
                                "bead_index": bead,
                                "type": "A" if types[chain, bead] == 1 else "B",
                                "type_numeric": int(types[chain, bead]),
                                "provenance": "deterministic reconstruction via repo generate_per_chain",
                            }
                        )
            for lag in range(0, 21):
                per_chain = np.mean(psi[:, : N_CHAIN - lag] * psi[:, lag:], axis=1)
                theory = 1.0 if lag == 0 else kappa**2 * (2.0 * PI_SELECTED - 1.0) ** lag
                cov_rows.append(
                    {
                        "kappa": kappa,
                        "pi": PI_SELECTED,
                        "f_A_target": F_A_TARGET,
                        "seed": seed,
                        "lag": lag,
                        "covariance_seed_mean": float(np.mean(per_chain)),
                        "covariance_chain_sd_within_seed": float(np.std(per_chain, ddof=1)),
                        "covariance_theory": theory,
                        "n_chains": N_CHAINS,
                        "provenance": "deterministic reconstruction; valid within-chain pairs only",
                    }
                )
            within = np.mean(psi[:, :-1] * psi[:, 1:], axis=1)
            across = psi[:-1, -1] * psi[1:, 0]
            for kind, values in (("within-chain lag 1", within), ("across-chain boundary", across)):
                boundary_rows.append(
                    {
                        "kappa": kappa,
                        "pi": PI_SELECTED,
                        "seed": seed,
                        "comparison": kind,
                        "covariance_seed_mean": float(np.mean(values)),
                        "covariance_unit_sd_within_seed": float(np.std(values, ddof=1)),
                        "n_units": int(len(values)),
                        "theory": (
                            kappa**2 * (2.0 * PI_SELECTED - 1.0)
                            if kind.startswith("within")
                            else 0.0
                        ),
                        "provenance": "deterministic reconstruction via repo generate_per_chain",
                    }
                )

    save(pd.DataFrame(strip_rows), "Fig1_sequence_strips.csv")
    save(pd.DataFrame(cov_rows), "Fig1_sequence_covariance.csv")
    save(pd.DataFrame(boundary_rows), "Fig1_chain_boundary_check.csv")

    # Analytic finite-chain zero-wavevector sequence form factor.
    rows = []
    for pi in (0.70, 0.85, 0.95, 0.99):
        for kappa in np.linspace(0.0, 1.0, 101):
            ell = np.arange(1, N_CHAIN)
            lam = 2.0 * pi - 1.0
            p0 = 1.0 + 2.0 * kappa**2 * np.sum((1.0 - ell / N_CHAIN) * lam**ell)
            rows.append(
                {
                    "pi": pi,
                    "kappa": kappa,
                    "lambda": lam,
                    "P_N_q0_normalized": p0,
                    "chain_length": N_CHAIN,
                }
            )
    save(pd.DataFrame(rows), "Fig1_finite_chain_form_factor_q0.csv")


def pi_kappa_tables() -> None:
    df = read_corrected("fig2_pi_kappa")
    df = df.drop(columns=["E", "xi"], errors="ignore")
    df["archive_to_number_scale"] = df.apply(archive_to_number_scale, axis=1)
    df["C_archive_N"] = df["C"] * df["archive_to_number_scale"]
    df["archive_final_frames_averaged"] = 5
    df["first_mixed_bin"] = np.isclose(df["k_star"], 0.4283989982167899)
    save(df, "Fig2_pi_kappa_runs.csv")

    rows = []
    counts = []
    for (pi, kappa), sub in df.groupby(["pi", "kappa"], sort=True):
        vc = sub["k_star"].value_counts()
        mode_k = float(vc[vc == vc.max()].index.min())
        rows.append(
            {
                "pi": pi,
                "kappa": kappa,
                "C_archive_N_mean": sub["C_archive_N"].mean(),
                "C_archive_N_sem": sem(sub["C_archive_N"]),
                "n_seeds": len(sub),
                "k_star_mean": sub["k_star"].mean(),
                "k_star_sem": sem(sub["k_star"]),
                "k_star_mode": mode_k,
                "first_mixed_bin_fraction": sub["first_mixed_bin"].mean(),
                    "metric_note": (
                        "number-normalized archived S_AA-S_AB peak; at f_A=0.5 its "
                    "ensemble mean is A/B-exchange symmetric; processed archive averages final 5 frames"
                    ),
            }
        )
        for q, n in vc.sort_index().items():
            counts.append(
                {
                    "pi": pi,
                    "kappa": kappa,
                    "k_star": q,
                    "n_runs": int(n),
                    "fraction": n / len(sub),
                }
            )
    save(pd.DataFrame(rows), "Fig2_pi_kappa_summary.csv")
    save(pd.DataFrame(counts), "Fig2_kstar_counts.csv")


def make_modes(n2_max: int = 49) -> tuple[np.ndarray, np.ndarray]:
    nmax = int(math.ceil(math.sqrt(n2_max)))
    modes = []
    n2 = []
    for nx in range(-nmax, nmax + 1):
        for ny in range(-nmax, nmax + 1):
            for nz in range(-nmax, nmax + 1):
                s = nx * nx + ny * ny + nz * nz
                if 0 < s <= n2_max:
                    modes.append((nx, ny, nz))
                    n2.append(s)
    return np.asarray(modes, dtype=int), np.asarray(n2, dtype=int)


def cic_grid(pos: np.ndarray, weights: np.ndarray, L: float, G: int) -> np.ndarray:
    h = L / G
    p = (np.mod(pos, L)) / h
    i0 = np.floor(p).astype(np.int64) % G
    frac = p - np.floor(p)
    i1 = (i0 + 1) % G
    wx0, wy0, wz0 = 1.0 - frac[:, 0], 1.0 - frac[:, 1], 1.0 - frac[:, 2]
    wx1, wy1, wz1 = frac[:, 0], frac[:, 1], frac[:, 2]
    out = np.zeros((G, G, G), dtype=np.float64)
    w = weights.astype(np.float64)
    for ix, iy, iz, xa, ya, za in (
        (i0[:, 0], i0[:, 1], i0[:, 2], wx0, wy0, wz0),
        (i1[:, 0], i0[:, 1], i0[:, 2], wx1, wy0, wz0),
        (i0[:, 0], i1[:, 1], i0[:, 2], wx0, wy1, wz0),
        (i0[:, 0], i0[:, 1], i1[:, 2], wx0, wy0, wz1),
        (i1[:, 0], i1[:, 1], i0[:, 2], wx1, wy1, wz0),
        (i1[:, 0], i0[:, 1], i1[:, 2], wx1, wy0, wz1),
        (i0[:, 0], i1[:, 1], i1[:, 2], wx0, wy1, wz1),
        (i1[:, 0], i1[:, 1], i1[:, 2], wx1, wy1, wz1),
    ):
        np.add.at(out, (ix, iy, iz), w * xa * ya * za)
    return out


def reciprocal_amplitudes(
    pos: np.ndarray, modes: np.ndarray, types: np.ndarray, L: float, chunk: int = 192
) -> tuple[np.ndarray, np.ndarray]:
    """Exact particle sums for A and B at reciprocal-lattice wavevectors."""
    wrapped = np.mod(pos, L).astype(np.float64, copy=False)
    q = (2.0 * np.pi / L) * modes.astype(float)
    a_mask = types == 1
    out_a = np.empty(len(modes), dtype=np.complex128)
    out_b = np.empty(len(modes), dtype=np.complex128)
    for start in range(0, len(modes), chunk):
        stop = min(start + chunk, len(modes))
        phase = np.exp(-1j * (wrapped @ q[start:stop].T))
        out_a[start:stop] = np.sum(phase[a_mask], axis=0)
        out_b[start:stop] = np.sum(phase[~a_mask], axis=0)
    return out_a, out_b


def spectra_tables() -> None:
    modes, n2 = make_modes(49)
    shell_values = sorted(np.unique(n2))
    per_frame: list[dict] = []
    metadata: list[dict] = []

    for kappa in (0.0, 0.5, 1.0):
        path = TRAJ / f"k{kappa:.1f}_tw100000000_s1_trajectory.npz"
        arr = np.load(path)
        positions = arr["positions"][-25:]
        steps = arr["steps"][-25:]
        types = arr["types"].astype(int)
        L = float(arr["box_size"])
        n_tot = len(types)
        metadata.append(
            {
                "kappa": kappa,
                "pi": PI_SELECTED,
                "f_A_target": F_A_TARGET,
                "f_A_realized": float(types.mean()),
                "seed": 1,
                "n_chains": N_CHAINS,
                "chain_length": N_CHAIN,
                "n_beads": n_tot,
                "box_size": L,
                "grid_size_for_CIC": G,
                "prequench_preparation_steps": 100000000,
                "postquench_total_steps": 250000,
                "late_frames_used": len(positions),
                "first_frame_step": int(steps[0]),
                "last_frame_step": int(steps[-1]),
                "bond_k": 100.0,
                "T_quench_star": 0.7,
                "eps_AA": 1.0,
                "eps_BB": 1.0,
                "eps_AB": 0.1,
                "n_reciprocal_modes": len(modes),
                "n_shells": len(shell_values),
                "q_max": K_FUND * math.sqrt(max(shell_values)),
                "scope_note": "representative single-seed late post-quench spectrum; not aging evidence",
            }
        )

        mode_indices = modes % G
        for frame_index, (step, pos) in enumerate(zip(steps, positions)):
            amp_a, amp_b = reciprocal_amplitudes(pos, modes, types, L)
            saa = np.abs(amp_a) ** 2 / n_tot
            sbb = np.abs(amp_b) ** 2 / n_tot
            sab = np.real(amp_a * np.conj(amp_b)) / n_tot
            spsi = saa + sbb - 2.0 * sab
            snn = saa + sbb + 2.0 * sab

            grid_a = cic_grid(pos, types.astype(float), L, G)
            grid_b = cic_grid(pos, 1.0 - types.astype(float), L, G)
            fft_a = np.fft.fftn(grid_a)
            fft_b = np.fft.fftn(grid_b)
            ia, ib, ic = mode_indices.T
            camp_a = fft_a[ia, ib, ic]
            camp_b = fft_b[ia, ib, ic]
            csaa = np.abs(camp_a) ** 2 / n_tot
            csbb = np.abs(camp_b) ** 2 / n_tot
            csab = np.real(camp_a * np.conj(camp_b)) / n_tot
            cspsi = csaa + csbb - 2.0 * csab
            csnn = csaa + csbb + 2.0 * csab

            # Standard first-order-assignment amplitude window.  The directly
            # summed spectrum remains the reference; this is only a diagnostic.
            window = np.prod(np.sinc(modes / G) ** 2, axis=1)
            da = camp_a / window
            db = camp_b / window
            dsaa = np.abs(da) ** 2 / n_tot
            dsbb = np.abs(db) ** 2 / n_tot
            dsab = np.real(da * np.conj(db)) / n_tot
            dspsi = dsaa + dsbb - 2.0 * dsab

            for shell in shell_values:
                mask = n2 == shell
                qval = K_FUND * math.sqrt(shell)
                row = {
                    "kappa": kappa,
                    "frame_index_within_late_window": frame_index,
                    "postquench_step": int(step),
                    "n2": shell,
                    "q_sigma": qval,
                    "n_modes": int(mask.sum()),
                    "S_AA_N_direct": float(np.mean(saa[mask])),
                    "S_BB_N_direct": float(np.mean(sbb[mask])),
                    "S_AB_N_direct": float(np.mean(sab[mask])),
                    "S_psipsi_N_direct": float(np.mean(spsi[mask])),
                    "C_N_direct": float(0.5 * np.mean(spsi[mask])),
                    "S_nn_N_direct": float(np.mean(snn[mask])),
                    "S_AA_N_CIC": float(np.mean(csaa[mask])),
                    "S_BB_N_CIC": float(np.mean(csbb[mask])),
                    "S_AB_N_CIC": float(np.mean(csab[mask])),
                    "S_psipsi_N_CIC": float(np.mean(cspsi[mask])),
                    "C_N_CIC": float(0.5 * np.mean(cspsi[mask])),
                    "S_nn_N_CIC": float(np.mean(csnn[mask])),
                    "C_N_CIC_window_deconvolved": float(0.5 * np.mean(dspsi[mask])),
                }
                per_frame.append(row)

    frame_df = pd.DataFrame(per_frame)
    save(frame_df, "Fig3_representative_spectra_per_frame.csv")
    save(pd.DataFrame(metadata), "Fig3_representative_spectra_metadata.csv")

    numeric = [
        c
        for c in frame_df.columns
        if c.startswith("S_") or c.startswith("C_N")
    ]
    rows = []
    for (kappa, n2v, q, nmodes), sub in frame_df.groupby(
        ["kappa", "n2", "q_sigma", "n_modes"], sort=True
    ):
        row = {
            "kappa": kappa,
            "n2": n2v,
            "q_sigma": q,
            "n_modes": nmodes,
            "n_late_frames": len(sub),
        }
        for col in numeric:
            row[f"{col}_mean"] = sub[col].mean()
            row[f"{col}_sd_over_frames"] = sub[col].std(ddof=1)
        row["CIC_raw_bias_percent"] = 100.0 * (
            row["C_N_CIC_mean"] / row["C_N_direct_mean"] - 1.0
        )
        row["CIC_deconvolved_bias_percent"] = 100.0 * (
            row["C_N_CIC_window_deconvolved_mean"] / row["C_N_direct_mean"] - 1.0
        )
        rows.append(row)
    save(pd.DataFrame(rows), "Fig3_representative_spectra_summary.csv")


def composition_tables() -> None:
    df = pd.read_csv(INPUT / "composition_symmetric_summary.csv")
    df["archive_to_number_scale"] = SCALE_DEFAULT
    df["C_symmetric_N"] = 0.5 * df["Spsipsi_peak"] * SCALE_DEFAULT
    df["symmetric_reanalysis_final_frames_averaged"] = 25
    df["C_symmetric_N_comp_normalized"] = df["C_symmetric_N"] / (
        4.0 * df["f_A"] * (1.0 - df["f_A"])
    )
    save(df, "Fig4_composition_runs_symmetric.csv")

    rows = []
    for (f_a, kappa), sub in df.groupby(["f_A", "kappa"], sort=True):
        rows.append(
            {
                "f_A": f_a,
                "kappa": kappa,
                "pi": float(sub["pi"].iloc[0]),
                "C_symmetric_N_mean": sub["C_symmetric_N"].mean(),
                "C_symmetric_N_sem": sem(sub["C_symmetric_N"]),
                "C_comp_normalized_mean": sub["C_symmetric_N_comp_normalized"].mean(),
                "C_comp_normalized_sem": sem(sub["C_symmetric_N_comp_normalized"]),
                "k_star_mean": sub["Spsipsi_k_star"].mean(),
                "k_star_sem": sem(sub["Spsipsi_k_star"]),
                "n_seeds": len(sub),
            }
        )
    summary = pd.DataFrame(rows)
    save(summary, "Fig4_composition_summary.csv")

    sym_rows = []
    for kappa in sorted(summary["kappa"].unique()):
        sub = summary[summary["kappa"] == kappa].set_index("f_A")
        for f_a in (0.2, 0.3, 0.4):
            left = sub.loc[f_a]
            right = sub.loc[1.0 - f_a]
            ratio = left["C_symmetric_N_mean"] / right["C_symmetric_N_mean"]
            rel = math.sqrt(
                (left["C_symmetric_N_sem"] / left["C_symmetric_N_mean"]) ** 2
                + (right["C_symmetric_N_sem"] / right["C_symmetric_N_mean"]) ** 2
            )
            sym_rows.append(
                {
                    "kappa": kappa,
                    "f_A": f_a,
                    "paired_f_A": 1.0 - f_a,
                    "C_ratio_f_over_1minusf": ratio,
                    "ratio_propagated_sem": ratio * rel,
                }
            )
    save(pd.DataFrame(sym_rows), "Fig4_composition_symmetry_ratios.csv")


def interaction_tables() -> None:
    df = read_corrected("fig4_epsAB")
    df = df.drop(columns=["E", "xi"], errors="ignore")
    df["archive_to_number_scale"] = df.apply(archive_to_number_scale, axis=1)
    df["C_archive_N"] = df["C"] * df["archive_to_number_scale"]
    df["archive_final_frames_averaged"] = 5
    save(df, "Fig5_interaction_runs.csv")

    rows = []
    for (sequence, eps), sub in df.groupby(["sequence", "eps_AB"], sort=True):
        rows.append(
            {
                "sequence": sequence,
                "eps_AB": eps,
                "C_archive_N_mean": sub["C_archive_N"].mean(),
                "C_archive_N_sem": sem(sub["C_archive_N"]),
                "n_seeds": len(sub),
                "k_star_mode": float(sub["k_star"].mode().min()),
            }
        )
    save(pd.DataFrame(rows), "Fig5_interaction_summary.csv")

    corr = df[df["sequence"] == "correlated"].copy()
    low = corr[corr["eps_AB"] <= 0.2].pivot(
        index="seed", columns="eps_AB", values="C_archive_N"
    )
    fr_archive = friedmanchisquare(*[low[c].values for c in low.columns])
    permutation_p, permutation_mcse, permutation_exceedances = blocked_friedman_permutation(
        low.to_numpy(), n_draws=1_000_000, seed=20260803, chunk_size=20_000
    )
    a = corr[np.isclose(corr["eps_AB"], 0.1)].set_index("seed")["C_archive_N"]
    b = corr[np.isclose(corr["eps_AB"], 0.8)].set_index("seed")["C_archive_N"]
    endpoint_differences = (a - b).dropna()
    endpoint_differences = endpoint_differences[endpoint_differences != 0]
    endpoint_positive = int((endpoint_differences > 0).sum())
    endpoint_n = int(len(endpoint_differences))
    endpoint_sign_p = float(
        binomtest(endpoint_positive, endpoint_n, p=0.5, alternative="two-sided").pvalue
    )
    test_rows = [
        {
            "analysis": "available processed archived estimator",
            "metric": "number-normalized archived S_AA-S_AB peak",
            "comparison": "eps_AB <= 0.2, 8 repeated conditions",
            "test": "1,000,000 within-seed permutations; RNG seed 20260803",
            "statistic": float(fr_archive.statistic),
            "p_value": permutation_p,
            "monte_carlo_se": permutation_mcse,
            "n_seeds": 5,
            "permutation_exceedances": permutation_exceedances,
            "n_permutations": 1_000_000,
            "rng_seed": 20260803,
            "provenance": "recomputed from supplied per-run CSV; independent np.random.Generator.permuted relabelling within each seed block in fixed 20,000-draw chunks",
        },
        {
            "analysis": "available processed archived estimator",
            "metric": "number-normalized archived S_AA-S_AB peak",
            "comparison": "eps_AB <= 0.2, 8 repeated conditions",
            "test": "Friedman asymptotic",
            "statistic": float(fr_archive.statistic),
            "p_value": float(fr_archive.pvalue),
            "monte_carlo_se": np.nan,
            "n_seeds": 5,
            "provenance": "recomputed from supplied processed per-run CSV",
        },
        {
            "analysis": "prior raw-spectrum audit",
            "metric": "number-normalized symmetric C=0.5 S_psipsi(k*)",
            "comparison": "eps_AB <= 0.2, 8 repeated conditions",
            "test": "Friedman statistic with 1,000,000 within-seed permutations",
            "statistic": 7.133333333333,
            "p_value": 0.436949,
            "monte_carlo_se": 0.000496,
            "n_seeds": 5,
            "provenance": "retained prior audit result; full raw archive unavailable",
        },
        {
            "analysis": "prior raw-spectrum audit",
            "metric": "number-normalized symmetric C=0.5 S_psipsi(k*)",
            "comparison": "eps_AB <= 0.2, 8 repeated conditions",
            "test": "Friedman asymptotic",
            "statistic": 7.133333333333,
            "p_value": 0.41513,
            "monte_carlo_se": np.nan,
            "n_seeds": 5,
            "provenance": "retained prior audit result; full raw archive unavailable",
        },
        {
            "analysis": "available processed archived estimator",
            "metric": "number-normalized archived S_AA-S_AB peak",
            "comparison": "eps_AB=0.1 versus 0.8, paired seeds",
            "test": "exact two-sided sign test; 5/5 signs concordant",
            "statistic": endpoint_positive,
            "p_value": endpoint_sign_p,
            "monte_carlo_se": np.nan,
            "n_seeds": 5,
            "provenance": "recomputed from supplied processed per-run CSV",
        },
        {
            "analysis": "prior raw-spectrum audit",
            "metric": "number-normalized symmetric C=0.5 S_psipsi(k*)",
            "comparison": "eps_AB=0.1 versus 0.8, paired seeds",
            "test": "exact two-sided sign test; 5/5 signs concordant",
            "statistic": 5,
            "p_value": 0.0625,
            "monte_carlo_se": np.nan,
            "n_seeds": 5,
            "provenance": "retained prior audit result; full raw archive unavailable",
        },
    ]
    save(pd.DataFrame(test_rows), "Fig5_interaction_tests.csv")

    prior_endpoints = pd.DataFrame(
        [
            {
                "eps_AB": 0.1,
                "C_symmetric_N_mean": 63.8,
                "C_symmetric_N_sem": 5.3,
                "n_seeds": 5,
                "provenance": "retained prior raw-spectrum audit",
            },
            {
                "eps_AB": 0.8,
                "C_symmetric_N_mean": 31.7,
                "C_symmetric_N_sem": 2.5,
                "n_seeds": 5,
                "provenance": "retained prior raw-spectrum audit",
            },
        ]
    )
    save(prior_endpoints, "Fig5_interaction_symmetric_endpoints_audit.csv")


def predictor(kappa: float, pi: float, q: float, N: int = 40, b: float = 1.0) -> float:
    ell = np.arange(1, N)
    lam = 2.0 * pi - 1.0
    return float(
        1.0
        + 2.0
        * kappa**2
        * np.sum((1.0 - ell / N) * lam**ell * np.exp(-(q * b) ** 2 * ell / 6.0))
    )


def theory_audit_tables() -> None:
    runs = pd.read_csv(OUT / "Fig2_pi_kappa_runs.csv")
    runs["P_N_kstar_seed"] = [
        predictor(k, p, q) for k, p, q in zip(runs.kappa, runs.pi, runs.k_star)
    ]
    runs["P_N_q0_seed"] = [predictor(k, p, 0.0) for k, p in zip(runs.kappa, runs.pi)]
    save(
        runs[["run", "pi", "kappa", "seed", "k_star", "C_archive_N", "P_N_kstar_seed", "P_N_q0_seed"]],
        "Fig6_predictor_seed_values.csv",
    )
    response = (
        runs.groupby(["pi", "kappa"], as_index=False)
        .agg(
            C_archive_N_mean=("C_archive_N", "mean"),
            C_archive_N_sem=("C_archive_N", sem),
            P_N_kstar_mean_over_seeds=("P_N_kstar_seed", "mean"),
            P_N_q0=("P_N_q0_seed", "mean"),
            n_seeds=("seed", "size"),
        )
    )
    response["inverse_C_archive_N"] = 1.0 / response["C_archive_N_mean"]
    response["inverse_P_N_kstar"] = 1.0 / response["P_N_kstar_mean_over_seeds"]
    response["inverse_P_N_q0"] = 1.0 / response["P_N_q0"]
    response["inactive_parameter_boundary"] = (response.kappa == 0) | (response.pi == 0.5)
    response["primary_regression_included"] = ~response["inactive_parameter_boundary"]
    save(response, "Fig6_predictor_points.csv")

    primary = response[(response.kappa > 0) & (response.pi > 0.5)]
    r_k = linregress(primary["inverse_P_N_kstar"], primary["inverse_C_archive_N"])
    r_0 = linregress(primary["inverse_P_N_q0"], primary["inverse_C_archive_N"])
    base_runs = runs[(runs.kappa == 0) | (runs.pi == 0.5)]
    collapsed = pd.concat(
        [
            primary,
            pd.DataFrame(
                [
                    {
                        "pi": np.nan,
                        "kappa": 0.0,
                        "C_archive_N_mean": base_runs.C_archive_N.mean(),
                        "C_archive_N_sem": sem(base_runs.C_archive_N),
                        "P_N_kstar_mean_over_seeds": base_runs.P_N_kstar_seed.mean(),
                        "P_N_q0": 1.0,
                        "n_seeds": len(base_runs),
                        "inverse_C_archive_N": 1.0 / base_runs.C_archive_N.mean(),
                        "inverse_P_N_kstar": 1.0 / base_runs.P_N_kstar_seed.mean(),
                        "inverse_P_N_q0": 1.0,
                        "primary_regression_included": False,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    crk = linregress(collapsed.inverse_P_N_kstar, collapsed.inverse_C_archive_N)
    cr0 = linregress(collapsed.inverse_P_N_q0, collapsed.inverse_C_archive_N)
    regression = pd.DataFrame(
        [
            {
                "dataset": "available archived balanced estimator; active interior primary",
                "predictor": "finite q=k_star",
                "r_squared": r_k.rvalue**2,
                "slope": r_k.slope,
                "intercept": r_k.intercept,
                "n_group_means": len(primary),
                "provenance": "recomputed from processed per-run CSV",
            },
            {
                "dataset": "available archived balanced estimator; active interior primary",
                "predictor": "q=0",
                "r_squared": r_0.rvalue**2,
                "slope": r_0.slope,
                "intercept": r_0.intercept,
                "n_group_means": len(primary),
                "provenance": "recomputed from processed per-run CSV",
            },
            {
                "dataset": "collapsed inactive-boundary baseline sensitivity",
                "predictor": "finite q=k_star",
                "r_squared": crk.rvalue**2,
                "slope": crk.slope,
                "intercept": crk.intercept,
                "n_group_means": len(collapsed),
                "provenance": "one collapsed inactive boundary plus 35 active conditions",
            },
            {
                "dataset": "collapsed inactive-boundary baseline sensitivity",
                "predictor": "q=0",
                "r_squared": cr0.rvalue**2,
                "slope": cr0.slope,
                "intercept": cr0.intercept,
                "n_group_means": len(collapsed),
                "provenance": "one collapsed inactive boundary plus 35 active conditions",
            },
        ]
    )
    save(regression, "Fig6_predictor_regression_audit.csv")

    bsens = []
    for b in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
        bruns = runs.copy()
        bruns["P"] = [predictor(k, p, q, b=b) for k, p, q in zip(bruns.kappa, bruns.pi, bruns.k_star)]
        bg = bruns.groupby(["pi", "kappa"], as_index=False).agg(C=("C_archive_N", "mean"), P=("P", "mean"))
        bg = bg[(bg.kappa > 0) & (bg.pi > 0.5)]
        lr = linregress(1.0 / bg.P, 1.0 / bg.C)
        bsens.append({"b_over_sigma": b, "finite_q_r_squared_active_interior": lr.rvalue**2, "q0_r_squared_active_interior": r_0.rvalue**2})
    save(pd.DataFrame(bsens), "FigS7_predictor_b_sensitivity.csv")

    form_rows = []
    for pi in (0.50, 0.85, 0.99):
        for kappa in (0.3, 0.6, 1.0):
            for q in np.linspace(0.0, 3.0, 151):
                form_rows.append(
                    {
                        "pi": pi,
                        "kappa": kappa,
                        "q_sigma": q,
                        "P_N_q_normalized": predictor(kappa, pi, q),
                        "chain_length": N_CHAIN,
                        "statistical_segment_length_b": 1.0,
                    }
                )
    save(pd.DataFrame(form_rows), "Fig6_finite_chain_form_factor.csv")

    w_tilde_q0, quadrature_error = force_equivalent_kernel_q0()
    rho = N_TOTAL / L_BOX**3
    epsilon_aa = 1.0
    epsilon_ab = 0.1
    temperature_star = 0.7
    chi_bare_q0 = -rho * (epsilon_aa - epsilon_ab) * w_tilde_q0 / temperature_star
    ideal_inverse_q0 = 1.0
    bare_denominator = ideal_inverse_q0 - chi_bare_q0 / 2.0
    kernel_provenance = (
        "recomputed by scipy.integrate.quad from w(r) with r_m=2^(1/6)sigma and "
        "r_c=2.5sigma; force branches and parameters verified against repo/melt/integrator.py"
    )
    rpa = pd.DataFrame(
        [
            {
                "quantity": "baseline number density rho",
                "value": rho,
                "units": "sigma^-3",
                "provenance": "recomputed as (144 chains x 40 beads)/(22 sigma)^3 from corrected static-scan parameters",
            },
            {
                "quantity": "force-equivalent w_tilde(0)",
                "value": w_tilde_q0,
                "units": "sigma^3",
                "provenance": kernel_provenance,
            },
            {
                "quantity": "w_tilde(0) quadrature error estimate",
                "value": quadrature_error,
                "units": "sigma^3",
                "provenance": "sum of scipy.integrate.quad absolute-error estimates for the two continuous branches",
            },
            {
                "quantity": "bare chi_bare(0)",
                "value": chi_bare_q0,
                "units": "dimensionless",
                "provenance": "recomputed as -rho*(epsilon_AA-epsilon_AB)*w_tilde(0)/T_star with epsilon_AA=1, epsilon_AB=0.1 and T_star=0.7",
            },
            {
                "quantity": "ideal inverse term at q=0 for kappa=0",
                "value": ideal_inverse_q0,
                "units": "dimensionless",
                "provenance": "P_N(0)=1 for kappa=0",
            },
            {
                "quantity": "interaction subtraction chi_bare(0)/2",
                "value": chi_bare_q0 / 2.0,
                "units": "dimensionless",
                "provenance": "recomputed from chi_bare(0); symmetric composition projection of the matrix closure",
            },
            {
                "quantity": "bare composition inverse stiffness",
                "value": bare_denominator,
                "units": "dimensionless",
                "provenance": "formal q-to-zero limit of 1-chi_bare(0)/2; failure of the uncalibrated bare homogeneous closure, not a measured spinodal",
            },
        ]
    )
    save(rpa, "Fig6_rpa_q0_diagnostic.csv")

    audit = pd.DataFrame(
        [
            {
                "audit_item": "CIC self-term-aware amplitude change across pi-kappa scan",
                "minimum": 0.437,
                "maximum": 4.251,
                "units": "percent",
                "affected_peak_calls": 1,
                "denominator_runs": 240,
                "provenance": "retained prior raw-spectrum audit",
            },
            {
                "audit_item": "CIC amplitude change for representative pi=0.99 conditions",
                "minimum": 0.437,
                "maximum": 2.493,
                "units": "percent",
                "affected_peak_calls": np.nan,
                "denominator_runs": np.nan,
                "provenance": "retained prior raw-spectrum audit",
            },
            {
                "audit_item": "archived asymmetric versus revised symmetric peak rule",
                "minimum": np.nan,
                "maximum": np.nan,
                "units": "peak calls",
                "affected_peak_calls": 4,
                "denominator_runs": 240,
                "provenance": "retained prior raw-spectrum audit",
            },
            {
                "audit_item": "revised symmetric selections of first admitted mixed bin",
                "minimum": np.nan,
                "maximum": np.nan,
                "units": "runs",
                "affected_peak_calls": 98,
                "denominator_runs": 240,
                "provenance": "retained prior raw-spectrum audit",
            },
            {
                "audit_item": "reciprocal modes in first admitted mixed bin",
                "minimum": np.nan,
                "maximum": np.nan,
                "units": "modes",
                "affected_peak_calls": 26,
                "denominator_runs": np.nan,
                "provenance": "exact mode enumeration for L=22,G=56",
            },
            {
                "audit_item": "modes in that bin below nominal center cutoff",
                "minimum": np.nan,
                "maximum": np.nan,
                "units": "modes",
                "affected_peak_calls": 18,
                "denominator_runs": np.nan,
                "provenance": "exact mode enumeration for L=22,G=56",
            },
        ]
    )
    save(audit, "Fig6_measurement_audit.csv")

    bin_modes = pd.DataFrame(
        [
            {
                "n_squared": n2,
                "q_sigma": K_FUND * math.sqrt(n2),
                "mode_count": count,
                "in_first_admitted_radial_bin": True,
                "below_nominal_center_cutoff": K_FUND * math.sqrt(n2) < 1.5 * K_FUND,
                "radial_bin_center": 0.4283989982167899,
                "radial_bin_lower_edge": K_FUND,
                "radial_bin_upper_edge": 2.0 * K_FUND,
                "nominal_center_cutoff": 1.5 * K_FUND,
            }
            for n2, count in ((1, 6), (2, 12), (3, 8))
        ]
    )
    save(bin_modes, "Fig6_first_bin_mode_inventory.csv")


def condition_controls() -> None:
    specs = [
        ("fig1_baseline", "baseline"),
        ("fig1_dense", "dense"),
        ("fig1_soft", "soft"),
        ("fig1_short_chain", "short chain"),
    ]
    runs = []
    summary = []
    for stem, label in specs:
        df = read_corrected(stem)
        df = df.drop(columns=["E", "xi"], errors="ignore")
        df["condition"] = label
        df["archive_to_number_scale"] = df.apply(archive_to_number_scale, axis=1)
        df["C_archive_N"] = df["C"] * df["archive_to_number_scale"]
        df["archive_final_frames_averaged"] = 5
        runs.append(df)
        for kappa, sub in df.groupby("kappa"):
            vc = sub["k_star"].value_counts()
            summary.append(
                {
                    "condition": label,
                    "kappa": kappa,
                    "C_archive_N_mean": sub["C_archive_N"].mean(),
                    "C_archive_N_sem": sem(sub["C_archive_N"]),
                    "k_star_mode": float(vc[vc == vc.max()].index.min()),
                    "k_star_mean": sub["k_star"].mean(),
                    "Rg_mean": sub["Rg"].mean(),
                    "Rg_sem": sem(sub["Rg"]),
                    "T_inst_mean": sub["T_inst"].mean(),
                    "T_inst_sem": sem(sub["T_inst"]),
                    "n_seeds": len(sub),
                    "n_chains": int(sub["n_chains"].iloc[0]),
                    "chain_length": int(sub["N"].iloc[0]),
                    "box_size": float(sub["L"].iloc[0]),
                    "rho": float(sub["rho"].iloc[0]),
                    "eps_AA": float(sub["eps_AA"].iloc[0]),
                    "eps_BB": float(sub["eps_BB"].iloc[0]),
                    "eps_AB": float(sub["eps_AB"].iloc[0]),
                    "eps_core": 1.0,
                    "pi": float(sub["pi"].iloc[0]),
                    "f_A": float(sub["f_A"].iloc[0]),
                    "T_quench_star": float(sub["T_quench_star"].iloc[0]),
                    "T_quench_K": float(sub["T_quench_K"].iloc[0]),
                    "bond_k": 200.0,
                    "prequench_equilibration_steps": 30000,
                    "postquench_steps": 250000,
                    "processed_final_frames_averaged": 5,
                }
            )
    save(pd.concat(runs, ignore_index=True), "FigS6_condition_control_runs.csv")
    save(pd.DataFrame(summary), "FigS6_condition_control_summary.csv")


def normalization_table() -> None:
    rows = []
    for stem in (
        "fig1_baseline",
        "fig1_dense",
        "fig1_soft",
        "fig1_short_chain",
        "fig2_pi_kappa",
        "fig3_fA_kappa",
        "fig4_epsAB",
    ):
        df = read_corrected(stem)
        first = df.iloc[0]
        rows.append(
            {
                "scan": stem,
                "grid_size_G": G,
                "n_chains": int(first["n_chains"]),
                "chain_length": int(first["N"]),
                "n_total": int(first["n_chains"] * first["N"]),
                "box_size_L": float(first["L"]),
                "volume": float(first["L"] ** 3),
                "archive_to_number_factor_G6_over_NV": archive_to_number_scale(first),
                "formula": "S_N(q)=G^6 S_code(q)/(N_total L^3)",
            }
        )
    save(pd.DataFrame(rows), "normalization_factors.csv")


def stage_original_processed_tables() -> None:
    """Preserve the supplied processed inputs outside the display-source layer.

    These originals contain archived E and xi columns that are deliberately not
    carried into any figure-input run table.
    """
    for path in sorted(CORRECTED.glob("*.csv")):
        shutil.copy2(path, RAW_WARNING / path.name)


def clean_generated_outputs() -> None:
    """Remove only files generated by this script before a reproducible rebuild."""
    for pattern in ("*.csv", "*.json"):
        for path in OUT.glob(pattern):
            path.unlink()
    for path in RAW_WARNING.glob("*.csv"):
        path.unlink()


def manifest() -> None:
    rows = []
    for path in sorted(OUT.glob("*.csv")):
        df = pd.read_csv(path)
        rows.append(
            {
                "file": path.name,
                "rows": len(df),
                "columns": len(df.columns),
                "column_names": " | ".join(df.columns),
            }
        )
    save(pd.DataFrame(rows), "source_data_manifest.csv")


def main() -> None:
    required = [
        CORRECTED / "fig2_pi_kappa.csv",
        INPUT / "composition_symmetric_summary.csv",
        TRAJ / "k0.0_tw100000000_s1_trajectory.npz",
        TRAJ / "k0.5_tw100000000_s1_trajectory.npz",
        TRAJ / "k1.0_tw100000000_s1_trajectory.npz",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))

    clean_generated_outputs()
    sequence_tables()
    pi_kappa_tables()
    spectra_tables()
    composition_tables()
    interaction_tables()
    theory_audit_tables()
    condition_controls()
    normalization_table()
    stage_original_processed_tables()
    manifest()

    summary = {
        "input_processed_rows": 789,
        "positive_quantitative_rows": 576,
        "composition_diagnostic_rows": 168,
        "excluded_prequench_duration_rows": 45,
        "grid_size": G,
        "default_archive_to_number_factor": SCALE_DEFAULT,
        "representative_spectrum_frames_per_condition": 25,
        "representative_spectrum_modes": int(len(make_modes(49)[0])),
        "repository_commit": "23fa2d03b4199d1d4e0c196f477a9569d10481b8",
    }
    (OUT / "rebuild_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
