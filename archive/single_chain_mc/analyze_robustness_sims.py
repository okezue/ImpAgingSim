#!/usr/bin/env python3
"""
Analyze robustness simulation outputs and generate all figures.

This script processes outputs from run_robustness_sims.sh and generates:
- A2: Finite-size scaling figures
- A3: Parameter scan figures (kappa, pi)
- B2: Lag-grid sensitivity figures
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from robustness_analysis import (
    extract_timescales_with_bootstrap,
    compute_multi_threshold_tau,
    compute_tstar_smooth,
    bootstrap_ci,
    fit_power_law_bootstrap,
)

def load_all_robustness_data(robustness_dir: str = "output/robustness") -> dict:
    """Load all robustness simulation CSV files."""
    data = {
        "size_scaling": [],
        "kappa_scan": [],
        "pi_scan": [],
        "dense_grid": [],
    }

    # A2: Size scaling (N=20,30,40)
    for fpath in glob.glob(os.path.join(robustness_dir, "results_N*_eps*_*.csv")):
        fname = os.path.basename(fpath)
        # Parse N from filename
        parts = fname.replace(".csv", "").split("_")
        N = int(parts[1].replace("N", ""))
        df = pd.read_csv(fpath)
        df["N_file"] = N  # Add N from filename in case column is missing
        data["size_scaling"].append(df)

    # A3: Kappa scan
    for fpath in glob.glob(os.path.join(robustness_dir, "results_kappa*_eps*_correlated.csv")):
        fname = os.path.basename(fpath)
        parts = fname.replace(".csv", "").split("_")
        kappa = float(parts[1].replace("kappa", ""))
        df = pd.read_csv(fpath)
        df["kappa_file"] = kappa
        data["kappa_scan"].append(df)

    # A3: Pi scan
    for fpath in glob.glob(os.path.join(robustness_dir, "results_pi*_eps*_correlated.csv")):
        fname = os.path.basename(fpath)
        parts = fname.replace(".csv", "").split("_")
        pi = float(parts[1].replace("pi", ""))
        df = pd.read_csv(fpath)
        df["pi_file"] = pi
        data["pi_scan"].append(df)

    # B2: Dense grid
    for fpath in glob.glob(os.path.join(robustness_dir, "results_dense_eps*_*.csv")):
        df = pd.read_csv(fpath)
        df["grid"] = "dense"
        data["dense_grid"].append(df)

    return data


def plot_size_scaling_chi4star(size_data: list, outdir: str):
    """A2: chi4* vs tw for different N values."""
    if not size_data:
        print("No size scaling data found")
        return

    df = pd.concat(size_data, ignore_index=True)

    # Use N column if available, otherwise N_file
    if "N" not in df.columns:
        df["N"] = df["N_file"]

    ts_df = extract_timescales_with_bootstrap(df, n_bootstrap=1000)

    N_values = sorted(ts_df["N"].unique())
    epsilons = sorted(ts_df["epsilon"].unique())

    fig, axes = plt.subplots(2, len(epsilons), figsize=(5*len(epsilons), 8))
    if len(epsilons) == 1:
        axes = axes.reshape(2, 1)

    markers = ['o', 's', '^', 'D']
    colors_N = plt.cm.viridis(np.linspace(0.2, 0.8, len(N_values)))

    for col, eps in enumerate(epsilons):
        for row, ens in enumerate(["iid", "correlated"]):
            ax = axes[row, col]

            for i, N in enumerate(N_values):
                subset = ts_df[(ts_df["ensemble"] == ens) &
                              (ts_df["epsilon"] == eps) &
                              (ts_df["N"] == N)]
                subset = subset[subset["tw"] > 0].sort_values("tw")

                if len(subset) == 0:
                    continue

                tw = subset["tw"].values
                chi4 = subset["chi4_star"].values
                chi4_lo = subset["chi4_star_lo"].values
                chi4_hi = subset["chi4_star_hi"].values

                ax.errorbar(tw, chi4, yerr=[chi4 - chi4_lo, chi4_hi - chi4],
                           marker=markers[i % len(markers)], linestyle='-',
                           color=colors_N[i], label=f"N={N}", markersize=5, capsize=2)

            ax.set_xscale("log")
            ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=10)
            ax.set_ylabel(r"$\chi_4^*$", fontsize=10)
            ax.set_title(f"{ens}, $\\epsilon$={eps}", fontsize=11)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_size_scaling_chi4star.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_size_scaling_chi4star.pdf"))
    plt.close()
    print("  Saved: Fig_size_scaling_chi4star")


def plot_size_scaling_separation(size_data: list, outdir: str):
    """A2: Ensemble separation vs tw for different N."""
    if not size_data:
        return

    df = pd.concat(size_data, ignore_index=True)
    if "N" not in df.columns:
        df["N"] = df["N_file"]

    ts_df = extract_timescales_with_bootstrap(df, n_bootstrap=1000)

    N_values = sorted(ts_df["N"].unique())
    epsilons = [e for e in sorted(ts_df["epsilon"].unique()) if e > 0]

    fig, axes = plt.subplots(1, len(epsilons), figsize=(6*len(epsilons), 5))
    if len(epsilons) == 1:
        axes = [axes]

    colors_N = plt.cm.Set1(np.linspace(0, 0.8, len(N_values)))

    for i, eps in enumerate(epsilons):
        ax = axes[i]

        for j, N in enumerate(N_values):
            # Get iid and correlated data
            iid = ts_df[(ts_df["ensemble"] == "iid") &
                       (ts_df["epsilon"] == eps) &
                       (ts_df["N"] == N)]
            corr = ts_df[(ts_df["ensemble"] == "correlated") &
                        (ts_df["epsilon"] == eps) &
                        (ts_df["N"] == N)]

            if len(iid) == 0 or len(corr) == 0:
                continue

            merged = pd.merge(iid, corr, on="tw", suffixes=("_iid", "_corr"))
            merged = merged[merged["tw"] > 0].sort_values("tw")

            if len(merged) == 0:
                continue

            tw = merged["tw"].values
            delta = (merged["chi4_star_iid"] - merged["chi4_star_corr"]) / merged["chi4_star_corr"] * 100

            ax.plot(tw, delta, marker='o', linestyle='-', color=colors_N[j],
                   label=f"N={N}", markersize=6)

        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.set_xscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=11)
        ax.set_ylabel(r"$\Delta\chi_4^* / \chi_{4,corr}^*$ (%)", fontsize=11)
        ax.set_title(f"$\\epsilon$={eps}", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_size_scaling_rel_separation.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_size_scaling_rel_separation.pdf"))
    plt.close()
    print("  Saved: Fig_size_scaling_rel_separation")


def plot_kappa_scan(kappa_data: list, outdir: str, tw_target: float = 10000):
    """A3: chi4* vs kappa."""
    if not kappa_data:
        print("No kappa scan data found")
        return

    df = pd.concat(kappa_data, ignore_index=True)
    if "kappa" not in df.columns or df["kappa"].nunique() <= 1:
        df["kappa"] = df["kappa_file"]

    ts_df = extract_timescales_with_bootstrap(df, n_bootstrap=1000)
    ts_df = ts_df[ts_df["tw"] == tw_target]

    if len(ts_df) == 0:
        print(f"No data for tw={tw_target}")
        return

    epsilons = sorted(ts_df["epsilon"].unique())

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    colors = plt.cm.Set1(np.linspace(0, 1, len(epsilons)))

    for i, eps in enumerate(epsilons):
        subset = ts_df[ts_df["epsilon"] == eps].sort_values("kappa")
        if len(subset) == 0:
            continue

        kappa = subset["kappa"].values
        chi4 = subset["chi4_star"].values
        chi4_lo = subset["chi4_star_lo"].values
        chi4_hi = subset["chi4_star_hi"].values

        ax.errorbar(kappa, chi4, yerr=[chi4 - chi4_lo, chi4_hi - chi4],
                   marker='o', linestyle='-', color=colors[i],
                   label=f"$\\epsilon$={eps}", markersize=8, capsize=3)

    ax.set_xlabel(r"$\kappa$", fontsize=12)
    ax.set_ylabel(r"$\chi_4^*$", fontsize=12)
    ax.set_title(f"$\\chi_4^*$ vs $\\kappa$ at $t_w$={int(tw_target)}", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"Fig_kappa_scan_chi4star_tw{int(tw_target)}.png"), dpi=200)
    plt.savefig(os.path.join(outdir, f"Fig_kappa_scan_chi4star_tw{int(tw_target)}.pdf"))
    plt.close()
    print(f"  Saved: Fig_kappa_scan_chi4star_tw{int(tw_target)}")


def plot_pi_scan(pi_data: list, outdir: str, tw_target: float = 10000):
    """A3: chi4* vs pi."""
    if not pi_data:
        print("No pi scan data found")
        return

    df = pd.concat(pi_data, ignore_index=True)
    if "pi" not in df.columns or df["pi"].nunique() <= 1:
        df["pi"] = df["pi_file"]

    ts_df = extract_timescales_with_bootstrap(df, n_bootstrap=1000)
    ts_df = ts_df[ts_df["tw"] == tw_target]

    if len(ts_df) == 0:
        print(f"No data for tw={tw_target}")
        return

    epsilons = sorted(ts_df["epsilon"].unique())

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    colors = plt.cm.Set1(np.linspace(0, 1, len(epsilons)))

    for i, eps in enumerate(epsilons):
        subset = ts_df[ts_df["epsilon"] == eps].sort_values("pi")
        if len(subset) == 0:
            continue

        pi = subset["pi"].values
        chi4 = subset["chi4_star"].values
        chi4_lo = subset["chi4_star_lo"].values
        chi4_hi = subset["chi4_star_hi"].values

        ax.errorbar(pi, chi4, yerr=[chi4 - chi4_lo, chi4_hi - chi4],
                   marker='s', linestyle='-', color=colors[i],
                   label=f"$\\epsilon$={eps}", markersize=8, capsize=3)

    ax.set_xlabel(r"$\pi$", fontsize=12)
    ax.set_ylabel(r"$\chi_4^*$", fontsize=12)
    ax.set_title(f"$\\chi_4^*$ vs $\\pi$ at $t_w$={int(tw_target)}", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"Fig_pi_scan_chi4star_tw{int(tw_target)}.png"), dpi=200)
    plt.savefig(os.path.join(outdir, f"Fig_pi_scan_chi4star_tw{int(tw_target)}.pdf"))
    plt.close()
    print(f"  Saved: Fig_pi_scan_chi4star_tw{int(tw_target)}")


def plot_grid_sensitivity(dense_data: list, baseline_csv: str, outdir: str):
    """B2: Compare coarse vs dense lag grid."""
    if not dense_data:
        print("No dense grid data found")
        return

    # Load baseline (coarse grid) data
    if not os.path.exists(baseline_csv):
        print(f"Baseline CSV not found: {baseline_csv}")
        return

    baseline = pd.read_csv(baseline_csv)
    dense = pd.concat(dense_data, ignore_index=True)

    # Extract timescales for both
    ts_coarse = extract_timescales_with_bootstrap(baseline, n_bootstrap=500)
    ts_coarse["grid"] = "coarse (28)"

    ts_dense = extract_timescales_with_bootstrap(dense, n_bootstrap=500)
    ts_dense["grid"] = "dense (56)"

    # Combine
    ts_all = pd.concat([ts_coarse, ts_dense], ignore_index=True)

    # Filter to common conditions
    common_tw = [1000, 3000, 10000]
    ts_all = ts_all[ts_all["tw"].isin(common_tw)]

    epsilons = sorted(ts_all["epsilon"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel A: t*_coarse vs t*_dense scatter
    ax = axes[0]

    coarse_vals = ts_all[ts_all["grid"] == "coarse (28)"].set_index(["ensemble", "epsilon", "tw"])
    dense_vals = ts_all[ts_all["grid"] == "dense (56)"].set_index(["ensemble", "epsilon", "tw"])

    common_idx = coarse_vals.index.intersection(dense_vals.index)

    if len(common_idx) > 0:
        t_coarse = coarse_vals.loc[common_idx, "t_star_discrete"].values
        t_dense = dense_vals.loc[common_idx, "t_star_discrete"].values

        ax.scatter(t_coarse, t_dense, alpha=0.7, s=50)
        lims = [min(t_coarse.min(), t_dense.min()), max(t_coarse.max(), t_dense.max())]
        ax.plot(lims, lims, 'k--', alpha=0.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t^*_{coarse}$ (discrete)", fontsize=11)
        ax.set_ylabel(r"$t^*_{dense}$ (discrete)", fontsize=11)
        ax.set_title("Discrete argmax", fontsize=12)
        ax.grid(True, alpha=0.3)

    # Panel B: same for smoothed
    ax = axes[1]
    if len(common_idx) > 0:
        t_coarse_s = coarse_vals.loc[common_idx, "t_star_smooth"].values
        t_dense_s = dense_vals.loc[common_idx, "t_star_smooth"].values

        ax.scatter(t_coarse_s, t_dense_s, alpha=0.7, s=50, color='green')
        lims = [min(t_coarse_s.min(), t_dense_s.min()), max(t_coarse_s.max(), t_dense_s.max())]
        ax.plot(lims, lims, 'k--', alpha=0.5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t^*_{coarse}$ (smooth)", fontsize=11)
        ax.set_ylabel(r"$t^*_{dense}$ (smooth)", fontsize=11)
        ax.set_title("Smoothed (quadratic)", fontsize=12)
        ax.grid(True, alpha=0.3)

    # Panel C: histogram of log ratio
    ax = axes[2]
    if len(common_idx) > 0:
        log_ratio_disc = np.log10(t_dense / t_coarse)
        log_ratio_smooth = np.log10(t_dense_s / t_coarse_s)

        ax.hist(log_ratio_disc, bins=15, alpha=0.5, label='Discrete', edgecolor='black')
        ax.hist(log_ratio_smooth, bins=15, alpha=0.5, label='Smooth', edgecolor='black')
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.7)
        ax.set_xlabel(r"$\log_{10}(t^*_{dense} / t^*_{coarse})$", fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.set_title("Grid sensitivity", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_tstar_grid_sensitivity.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_tstar_grid_sensitivity.pdf"))
    plt.close()
    print("  Saved: Fig_tstar_grid_sensitivity")


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--robustness_dir", default="output/robustness",
                  help="Directory containing robustness simulation outputs")
    p.add_argument("--baseline_csv", default="output/results.csv",
                  help="Baseline results CSV")
    p.add_argument("--outdir", default="analysis/robustness_full",
                  help="Output directory for figures")

    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading robustness simulation data...")
    data = load_all_robustness_data(args.robustness_dir)

    print(f"  Size scaling files: {len(data['size_scaling'])}")
    print(f"  Kappa scan files: {len(data['kappa_scan'])}")
    print(f"  Pi scan files: {len(data['pi_scan'])}")
    print(f"  Dense grid files: {len(data['dense_grid'])}")
    print()

    print("Generating figures...")

    # A2: Size scaling
    if data["size_scaling"]:
        plot_size_scaling_chi4star(data["size_scaling"], args.outdir)
        plot_size_scaling_separation(data["size_scaling"], args.outdir)

    # A3: Parameter scans
    if data["kappa_scan"]:
        plot_kappa_scan(data["kappa_scan"], args.outdir, tw_target=10000)
        plot_kappa_scan(data["kappa_scan"], args.outdir, tw_target=1000)

    if data["pi_scan"]:
        plot_pi_scan(data["pi_scan"], args.outdir, tw_target=10000)
        plot_pi_scan(data["pi_scan"], args.outdir, tw_target=1000)

    # B2: Grid sensitivity
    if data["dense_grid"]:
        plot_grid_sensitivity(data["dense_grid"], args.baseline_csv, args.outdir)

    print()
    print("=" * 50)
    print("Robustness figure generation complete!")
    print("=" * 50)
    print(f"Figures saved to: {args.outdir}")


if __name__ == "__main__":
    main()
