#!/usr/bin/env python3
"""
Enhanced robustness analysis for IMP aging simulation.

Implements:
A1) Multi-threshold tau_c with right-censoring
A2) Finite-size analysis (group by N)
A3) Parameter scan analysis (kappa, pi)
B1) Smoothed t* estimator (quadratic interpolation in log-space)
B2) Lag-grid sensitivity comparison
B3) Bootstrap confidence intervals
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import linregress
from typing import Dict, List, Tuple, Optional
import warnings

# -----------------------------------------------------------------------------
# B1: Smoothed t* estimator
# -----------------------------------------------------------------------------

def compute_tstar_smooth(lags: np.ndarray, chi4_vals: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Compute smoothed t* using quadratic interpolation in log-space.

    Returns: (t_star_discrete, t_star_smooth, chi4_star_discrete, chi4_star_smooth)
    """
    if len(lags) < 3 or len(chi4_vals) < 3:
        return np.nan, np.nan, np.nan, np.nan

    # Filter positive lags
    mask = lags > 0
    lags = lags[mask]
    chi4_vals = chi4_vals[mask]

    if len(lags) < 3:
        return np.nan, np.nan, np.nan, np.nan

    # Discrete maximum
    k = np.argmax(chi4_vals)
    t_star_discrete = lags[k]
    chi4_star_discrete = chi4_vals[k]

    # Default smooth values to discrete
    t_star_smooth = t_star_discrete
    chi4_star_smooth = chi4_star_discrete

    # Try quadratic interpolation if not at boundary
    if 1 <= k <= len(lags) - 2:
        # Points in log-space
        x = np.log(lags[[k-1, k, k+1]])
        y = chi4_vals[[k-1, k, k+1]]

        # Fit quadratic: y = ax^2 + bx + c
        # Using Vandermonde matrix
        V = np.vstack([x**2, x, np.ones_like(x)]).T
        try:
            coeffs = np.linalg.lstsq(V, y, rcond=None)[0]
            a, b, c = coeffs

            # Vertex at x* = -b/(2a)
            if a < 0:  # Proper maximum
                x_star = -b / (2 * a)
                if x[0] <= x_star <= x[2]:
                    t_star_smooth = np.exp(x_star)
                    chi4_star_smooth = a * x_star**2 + b * x_star + c
        except:
            pass

    return t_star_discrete, t_star_smooth, chi4_star_discrete, chi4_star_smooth


# -----------------------------------------------------------------------------
# A1: Multi-threshold tau_c with right-censoring
# -----------------------------------------------------------------------------

def compute_tau_c(lags: np.ndarray, Q_norm: np.ndarray, c: float, t_max: float) -> Tuple[float, bool]:
    """
    Compute tau_c (first crossing time through threshold c).

    Returns: (tau_c, is_censored)
        - tau_c: crossing time or t_max if censored
        - is_censored: True if Q_norm never crossed c
    """
    if len(lags) < 2 or len(Q_norm) < 2:
        return np.nan, True

    # Find first crossing below c
    below = Q_norm < c
    if not np.any(below):
        return t_max, True  # Censored

    # Linear interpolation for more precise crossing
    idx = np.argmax(below)
    if idx == 0:
        return lags[0], False

    # Interpolate between idx-1 and idx
    x0, x1 = lags[idx-1], lags[idx]
    y0, y1 = Q_norm[idx-1], Q_norm[idx]

    if y0 == y1:
        return x1, False

    # Linear interpolation: find x where y = c
    tau_c = x0 + (c - y0) * (x1 - x0) / (y1 - y0)
    return tau_c, False


def compute_multi_threshold_tau(df: pd.DataFrame, thresholds: List[float] = [0.8, 0.6, np.exp(-1)]) -> pd.DataFrame:
    """
    Compute tau_c for multiple thresholds across all conditions.
    """
    results = []

    group_cols = ["ensemble", "epsilon"]
    if "N" in df.columns:
        group_cols.append("N")
    if "kappa" in df.columns:
        group_cols.append("kappa")
    if "pi" in df.columns:
        group_cols.append("pi")

    for keys, group in df.groupby(group_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_dict = dict(zip(group_cols, keys))

        # Aggregate over disorder
        agg = group.groupby(["tw", "lag"]).agg(
            Q_mean=("Q_mean", "mean"),
            chi4=("chi4", "mean"),
        ).reset_index()

        t_max = agg["lag"].max()

        for tw in agg["tw"].unique():
            tw_data = agg[agg["tw"] == tw].sort_values("lag")
            lags = tw_data["lag"].values
            Q_mean = tw_data["Q_mean"].values

            # Normalize Q
            Q0 = Q_mean[0] if len(Q_mean) > 0 and Q_mean[0] > 0 else 1.0
            Q_norm = Q_mean / Q0

            row = {**key_dict, "tw": tw}

            for c in thresholds:
                tau, censored = compute_tau_c(lags, Q_norm, c, t_max)
                c_label = f"{c:.2f}".replace(".", "p")
                row[f"tau_{c_label}"] = tau
                row[f"censored_{c_label}"] = censored

            results.append(row)

    return pd.DataFrame(results)


# -----------------------------------------------------------------------------
# B3: Bootstrap confidence intervals
# -----------------------------------------------------------------------------

def bootstrap_ci(values: np.ndarray, n_bootstrap: int = 2000, ci: float = 0.95) -> Tuple[float, float, float]:
    """
    Compute bootstrap confidence interval.

    Returns: (median, ci_low, ci_high)
    """
    if len(values) < 2:
        return np.nanmedian(values), np.nan, np.nan

    rng = np.random.default_rng(42)
    boot_medians = []

    for _ in range(n_bootstrap):
        sample = rng.choice(values, size=len(values), replace=True)
        boot_medians.append(np.nanmedian(sample))

    boot_medians = np.array(boot_medians)
    alpha = (1 - ci) / 2
    ci_low = np.nanpercentile(boot_medians, 100 * alpha)
    ci_high = np.nanpercentile(boot_medians, 100 * (1 - alpha))

    return np.nanmedian(values), ci_low, ci_high


def extract_timescales_with_bootstrap(df: pd.DataFrame, n_bootstrap: int = 2000) -> pd.DataFrame:
    """
    Extract t*, chi4* with bootstrap CIs over disorder realizations.
    """
    results = []

    group_cols = ["ensemble", "epsilon"]
    if "N" in df.columns:
        group_cols.append("N")
    # Include kappa and pi if they have multiple unique values (parameter scans)
    if "kappa" in df.columns and df["kappa"].nunique() > 1:
        group_cols.append("kappa")
    if "pi" in df.columns and df["pi"].nunique() > 1:
        group_cols.append("pi")

    for keys, group in df.groupby(group_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_dict = dict(zip(group_cols, keys))

        tws = group["tw"].unique()

        for tw in tws:
            tw_group = group[group["tw"] == tw]
            disorder_indices = tw_group["disorder_idx"].unique()

            # Collect t* and chi4* per disorder
            t_star_discrete_list = []
            t_star_smooth_list = []
            chi4_star_list = []

            for d_idx in disorder_indices:
                d_data = tw_group[tw_group["disorder_idx"] == d_idx].sort_values("lag")
                lags = d_data["lag"].values
                chi4_vals = d_data["chi4"].values

                t_disc, t_smooth, chi4_disc, chi4_smooth = compute_tstar_smooth(lags, chi4_vals)

                t_star_discrete_list.append(t_disc)
                t_star_smooth_list.append(t_smooth)
                chi4_star_list.append(chi4_smooth)

            # Bootstrap CIs
            t_disc_arr = np.array(t_star_discrete_list)
            t_smooth_arr = np.array(t_star_smooth_list)
            chi4_arr = np.array(chi4_star_list)

            t_disc_med, t_disc_lo, t_disc_hi = bootstrap_ci(t_disc_arr, n_bootstrap)
            t_smooth_med, t_smooth_lo, t_smooth_hi = bootstrap_ci(t_smooth_arr, n_bootstrap)
            chi4_med, chi4_lo, chi4_hi = bootstrap_ci(chi4_arr, n_bootstrap)

            results.append({
                **key_dict,
                "tw": tw,
                "t_star_discrete": t_disc_med,
                "t_star_discrete_lo": t_disc_lo,
                "t_star_discrete_hi": t_disc_hi,
                "t_star_smooth": t_smooth_med,
                "t_star_smooth_lo": t_smooth_lo,
                "t_star_smooth_hi": t_smooth_hi,
                "chi4_star": chi4_med,
                "chi4_star_lo": chi4_lo,
                "chi4_star_hi": chi4_hi,
                "n_disorder": len(disorder_indices),
            })

    return pd.DataFrame(results)


# -----------------------------------------------------------------------------
# Power-law fitting with bootstrap
# -----------------------------------------------------------------------------

def fit_power_law_bootstrap(x: np.ndarray, y: np.ndarray, n_bootstrap: int = 2000) -> Dict:
    """
    Fit power law y = A * x^nu with bootstrap CI for the exponent.
    """
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if len(x) < 3:
        return {"A": np.nan, "nu": np.nan, "nu_lo": np.nan, "nu_hi": np.nan, "R2": np.nan}

    log_x = np.log10(x)
    log_y = np.log10(y)

    # Direct fit
    slope, intercept, r_value, _, _ = linregress(log_x, log_y)

    # Bootstrap
    rng = np.random.default_rng(42)
    boot_slopes = []

    for _ in range(n_bootstrap):
        idx = rng.choice(len(x), size=len(x), replace=True)
        # Check for unique x values
        if len(np.unique(log_x[idx])) < 2:
            continue
        try:
            s, _, _, _, _ = linregress(log_x[idx], log_y[idx])
            boot_slopes.append(s)
        except:
            continue

    boot_slopes = np.array(boot_slopes)
    nu_lo = np.nanpercentile(boot_slopes, 2.5)
    nu_hi = np.nanpercentile(boot_slopes, 97.5)

    return {
        "A": 10**intercept,
        "nu": slope,
        "nu_lo": nu_lo,
        "nu_hi": nu_hi,
        "R2": r_value**2,
    }


# -----------------------------------------------------------------------------
# Plotting functions
# -----------------------------------------------------------------------------

def plot_tstar_discrete_vs_smooth(ts_df: pd.DataFrame, outdir: str):
    """B1: Scatter of discrete vs smooth t*."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Scatter plot
    ax = axes[0]
    valid = ts_df.dropna(subset=["t_star_discrete", "t_star_smooth"])
    ax.scatter(valid["t_star_discrete"], valid["t_star_smooth"], alpha=0.6, s=40)

    # Add y=x line
    lims = [min(valid["t_star_discrete"].min(), valid["t_star_smooth"].min()),
            max(valid["t_star_discrete"].max(), valid["t_star_smooth"].max())]
    ax.plot(lims, lims, 'k--', alpha=0.5, label="y=x")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t^*_{discrete}$", fontsize=11)
    ax.set_ylabel(r"$t^*_{smooth}$", fontsize=11)
    ax.set_title("Discrete vs Smoothed Peak Location")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Histogram of log ratio
    ax = axes[1]
    log_ratio = np.log10(valid["t_star_smooth"] / valid["t_star_discrete"])
    ax.hist(log_ratio, bins=20, edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel(r"$\log_{10}(t^*_{smooth} / t^*_{discrete})$", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("Distribution of Smoothing Effect")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_tstar_discrete_vs_smooth_scatter.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_tstar_discrete_vs_smooth_scatter.pdf"))
    plt.close()


def plot_tstar_with_ci(ts_df: pd.DataFrame, outdir: str):
    """B3: t* vs tw with bootstrap confidence intervals."""
    epsilons = sorted(ts_df["epsilon"].unique())
    n_eps = len(epsilons)

    fig, axes = plt.subplots(1, n_eps, figsize=(5*n_eps, 4.5))
    if n_eps == 1:
        axes = [axes]

    for i, eps in enumerate(epsilons):
        ax = axes[i]
        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            subset = ts_df[(ts_df["ensemble"] == ens) & (ts_df["epsilon"] == eps)]
            subset = subset[subset["tw"] > 0].sort_values("tw")

            if len(subset) == 0:
                continue

            tw = subset["tw"].values
            t_star = subset["t_star_smooth"].values
            t_lo = subset["t_star_smooth_lo"].values
            t_hi = subset["t_star_smooth_hi"].values

            ax.errorbar(tw, t_star, yerr=[t_star - t_lo, t_hi - t_star],
                       color=color, marker=marker, linestyle='-', capsize=3,
                       label=ens, markersize=6)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=11)
        ax.set_ylabel(r"$t^*_{smooth}$", fontsize=11)
        ax.set_title(f"$\\epsilon$ = {eps}", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_tstar_vs_tw_with_CI.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_tstar_vs_tw_with_CI.pdf"))
    plt.close()


def plot_tau_c_thresholds(tau_df: pd.DataFrame, outdir: str):
    """A1: tau_c for multiple thresholds."""
    epsilons = sorted(tau_df["epsilon"].unique())
    ensembles = ["iid", "correlated"]

    thresholds = [0.8, 0.6, np.exp(-1)]
    threshold_labels = ["0.8", "0.6", r"$e^{-1}$"]
    colors = ["green", "orange", "purple"]

    fig, axes = plt.subplots(2, len(epsilons), figsize=(5*len(epsilons), 8))
    if len(epsilons) == 1:
        axes = axes.reshape(2, 1)

    for col, eps in enumerate(epsilons):
        for row, ens in enumerate(ensembles):
            ax = axes[row, col]
            subset = tau_df[(tau_df["ensemble"] == ens) & (tau_df["epsilon"] == eps)]
            subset = subset[subset["tw"] > 0].sort_values("tw")

            for c, label, color in zip(thresholds, threshold_labels, colors):
                c_label = f"{c:.2f}".replace(".", "p")
                col_name = f"tau_{c_label}"
                cens_col = f"censored_{c_label}"

                if col_name not in subset.columns:
                    continue

                tw = subset["tw"].values
                tau = subset[col_name].values
                censored = subset[cens_col].values if cens_col in subset.columns else np.zeros(len(tw), dtype=bool)

                # Plot uncensored points
                uncens_mask = ~censored
                if np.any(uncens_mask):
                    ax.plot(tw[uncens_mask], tau[uncens_mask], marker='o', linestyle='-',
                           color=color, label=f"c={label}", markersize=5)

                # Plot censored points with different marker
                if np.any(censored):
                    ax.scatter(tw[censored], tau[censored], marker='^', facecolors='none',
                              edgecolors=color, s=50, label=f"c={label} (censored)" if np.any(uncens_mask) else f"c={label}")

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=10)
            ax.set_ylabel(r"$\tau_c$", fontsize=10)
            ax.set_title(f"{ens}, $\\epsilon$={eps}", fontsize=11)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_tau_c_thresholds_vs_tw.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_tau_c_thresholds_vs_tw.pdf"))
    plt.close()


def plot_size_scaling(ts_df: pd.DataFrame, outdir: str):
    """A2: chi4* vs tw for different N values."""
    if "N" not in ts_df.columns:
        print("No N column found, skipping size scaling plot")
        return

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

                ax.plot(tw, chi4, marker=markers[i % len(markers)], linestyle='-',
                       color=colors_N[i], label=f"N={N}", markersize=5)

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


def plot_kappa_scan(ts_df: pd.DataFrame, outdir: str, tw_target: float = 10000):
    """A3: chi4* vs kappa for parameter scan."""
    if "kappa" not in ts_df.columns:
        print("No kappa column found, skipping kappa scan plot")
        return

    # Filter to target tw
    ts_df = ts_df[ts_df["tw"] == tw_target]

    if len(ts_df) == 0:
        print(f"No data for tw={tw_target}")
        return

    epsilons = sorted(ts_df["epsilon"].unique())

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    colors = plt.cm.Set1(np.linspace(0, 1, len(epsilons)))

    for i, eps in enumerate(epsilons):
        subset = ts_df[(ts_df["ensemble"] == "correlated") & (ts_df["epsilon"] == eps)]
        subset = subset.sort_values("kappa")

        if len(subset) == 0:
            continue

        kappa = subset["kappa"].values
        chi4 = subset["chi4_star"].values

        ax.plot(kappa, chi4, marker='o', linestyle='-', color=colors[i],
               label=f"$\\epsilon$={eps}", markersize=8)

    ax.set_xlabel(r"$\kappa$", fontsize=12)
    ax.set_ylabel(r"$\chi_4^*$", fontsize=12)
    ax.set_title(f"$\\chi_4^*$ vs $\\kappa$ at $t_w$={int(tw_target)}", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"Fig_kappa_scan_chi4star_tw{int(tw_target)}.png"), dpi=200)
    plt.savefig(os.path.join(outdir, f"Fig_kappa_scan_chi4star_tw{int(tw_target)}.pdf"))
    plt.close()


def plot_nu_bootstrap(ts_df: pd.DataFrame, outdir: str, n_bootstrap: int = 2000):
    """B3: nu (scaling exponent) with bootstrap CI."""
    epsilons = sorted(ts_df["epsilon"].unique())
    ensembles = ["iid", "correlated"]

    results = []

    for eps in epsilons:
        for ens in ensembles:
            subset = ts_df[(ts_df["ensemble"] == ens) & (ts_df["epsilon"] == eps)]
            subset = subset[subset["tw"] > 0].sort_values("tw")

            if len(subset) < 3:
                continue

            tw = subset["tw"].values
            t_star = subset["t_star_smooth"].values

            fit = fit_power_law_bootstrap(tw, t_star, n_bootstrap)

            results.append({
                "epsilon": eps,
                "ensemble": ens,
                **fit
            })

    fit_df = pd.DataFrame(results)

    if len(fit_df) == 0:
        print("No fit results to plot")
        return fit_df

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: nu vs epsilon
    ax = axes[0]
    for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
        subset = fit_df[fit_df["ensemble"] == ens]
        if len(subset) == 0:
            continue

        ax.errorbar(subset["epsilon"], subset["nu"],
                   yerr=[subset["nu"] - subset["nu_lo"], subset["nu_hi"] - subset["nu"]],
                   color=color, marker=marker, linestyle='-', capsize=4, markersize=8,
                   label=ens)

    ax.set_xlabel(r"$\epsilon$ (disorder strength)", fontsize=11)
    ax.set_ylabel(r"$\nu$ (exponent in $t^* \sim t_w^\nu$)", fontsize=11)
    ax.set_title("Scaling exponent with 95% CI", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: Delta nu
    ax = axes[1]
    iid_fits = fit_df[fit_df["ensemble"] == "iid"].set_index("epsilon")
    corr_fits = fit_df[fit_df["ensemble"] == "correlated"].set_index("epsilon")

    common_eps = [e for e in epsilons if e in iid_fits.index and e in corr_fits.index]
    delta_nu = []
    delta_nu_err = []

    for eps in common_eps:
        nu_iid = iid_fits.loc[eps, "nu"]
        nu_corr = corr_fits.loc[eps, "nu"]

        # Propagate errors (approximate)
        err_iid = (iid_fits.loc[eps, "nu_hi"] - iid_fits.loc[eps, "nu_lo"]) / 4
        err_corr = (corr_fits.loc[eps, "nu_hi"] - corr_fits.loc[eps, "nu_lo"]) / 4

        delta_nu.append(nu_iid - nu_corr)
        delta_nu_err.append(np.sqrt(err_iid**2 + err_corr**2))

    ax.errorbar(common_eps, delta_nu, yerr=delta_nu_err,
               color='purple', marker='D', linestyle='-', capsize=4, markersize=8)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel(r"$\epsilon$ (disorder strength)", fontsize=11)
    ax.set_ylabel(r"$\Delta\nu = \nu_{iid} - \nu_{corr}$", fontsize=11)
    ax.set_title("Exponent difference with 95% CI", fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_nu_star_bootstrap.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_nu_star_bootstrap.pdf"))
    plt.close()

    return fit_df


def plot_chi4_peak_fit_example(df: pd.DataFrame, outdir: str,
                               ensemble: str = "iid", epsilon: float = 3.0,
                               tw: float = 3000):
    """B1: Example of quadratic fit to chi4 peak."""
    # Get aggregated data for this condition
    subset = df[(df["ensemble"] == ensemble) & (df["epsilon"] == epsilon) & (df["tw"] == tw)]

    if len(subset) == 0:
        print(f"No data for {ensemble}, eps={epsilon}, tw={tw}")
        return

    agg = subset.groupby("lag").agg(chi4=("chi4", "mean")).reset_index()
    agg = agg.sort_values("lag")

    lags = agg["lag"].values
    chi4 = agg["chi4"].values

    # Filter positive lags
    mask = lags > 0
    lags = lags[mask]
    chi4 = chi4[mask]

    # Find discrete max
    k = np.argmax(chi4)
    t_star_discrete = lags[k]

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    # Plot data
    ax.plot(lags, chi4, 'ko-', label='Data', markersize=6)

    # Quadratic fit if not at boundary
    if 1 <= k <= len(lags) - 2:
        x_fit = np.log(lags[[k-1, k, k+1]])
        y_fit = chi4[[k-1, k, k+1]]

        V = np.vstack([x_fit**2, x_fit, np.ones_like(x_fit)]).T
        coeffs = np.linalg.lstsq(V, y_fit, rcond=None)[0]
        a, b, c = coeffs

        # Plot fit curve
        x_dense = np.linspace(x_fit[0], x_fit[-1], 100)
        y_dense = a * x_dense**2 + b * x_dense + c
        ax.plot(np.exp(x_dense), y_dense, 'r-', linewidth=2, label='Quadratic fit')

        # Highlight fit points
        ax.scatter(np.exp(x_fit), y_fit, color='red', s=100, zorder=5, marker='s')

        # Vertex
        if a < 0:
            x_star = -b / (2*a)
            t_star_smooth = np.exp(x_star)
            chi4_star_smooth = a * x_star**2 + b * x_star + c

            ax.axvline(x=t_star_discrete, color='blue', linestyle='--', alpha=0.7,
                      label=f'$t^*_{{discrete}}$ = {t_star_discrete:.0f}')
            ax.axvline(x=t_star_smooth, color='green', linestyle='-', alpha=0.7,
                      label=f'$t^*_{{smooth}}$ = {t_star_smooth:.1f}')

    ax.set_xscale("log")
    ax.set_xlabel("lag (MC sweeps)", fontsize=11)
    ax.set_ylabel(r"$\chi_4$", fontsize=11)
    ax.set_title(f"Peak fit example: {ensemble}, $\\epsilon$={epsilon}, $t_w$={tw}", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_chi4_peak_fit_example.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_chi4_peak_fit_example.pdf"))
    plt.close()


# -----------------------------------------------------------------------------
# Main analysis function
# -----------------------------------------------------------------------------

def run_robustness_analysis(csv_path: str, outdir: str):
    """Run complete robustness analysis on simulation data."""

    os.makedirs(outdir, exist_ok=True)

    print(f"Loading data from {csv_path}")
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows")

    # 1. Extract timescales with bootstrap CIs
    print("\nExtracting timescales with bootstrap CIs...")
    ts_df = extract_timescales_with_bootstrap(df, n_bootstrap=2000)
    ts_df.to_csv(os.path.join(outdir, "timescales_bootstrap.csv"), index=False)
    print(f"Saved timescales to {outdir}/timescales_bootstrap.csv")

    # 2. Multi-threshold tau_c
    print("\nComputing multi-threshold tau_c...")
    tau_df = compute_multi_threshold_tau(df)
    tau_df.to_csv(os.path.join(outdir, "tau_thresholds.csv"), index=False)
    print(f"Saved tau thresholds to {outdir}/tau_thresholds.csv")

    # 3. Generate plots
    print("\nGenerating plots...")

    # B1: Discrete vs smooth t*
    plot_tstar_discrete_vs_smooth(ts_df, outdir)
    print("  - Fig_tstar_discrete_vs_smooth_scatter")

    # B1: Peak fit example
    plot_chi4_peak_fit_example(df, outdir, ensemble="iid", epsilon=3.0, tw=3000)
    print("  - Fig_chi4_peak_fit_example")

    # B3: t* with CI
    plot_tstar_with_ci(ts_df, outdir)
    print("  - Fig_tstar_vs_tw_with_CI")

    # A1: Multi-threshold tau
    plot_tau_c_thresholds(tau_df, outdir)
    print("  - Fig_tau_c_thresholds_vs_tw")

    # A2: Size scaling (if N column exists)
    if "N" in ts_df.columns and ts_df["N"].nunique() > 1:
        plot_size_scaling(ts_df, outdir)
        print("  - Fig_size_scaling_chi4star")

    # A3: Kappa scan (if kappa column varies)
    if "kappa" in ts_df.columns and ts_df["kappa"].nunique() > 1:
        plot_kappa_scan(ts_df, outdir)
        print("  - Fig_kappa_scan_chi4star")

    # B3: Bootstrap nu
    fit_df = plot_nu_bootstrap(ts_df, outdir)
    if len(fit_df) > 0:
        fit_df.to_csv(os.path.join(outdir, "powerlaw_fits_bootstrap.csv"), index=False)
        print("  - Fig_nu_star_bootstrap")

    print("\n" + "="*60)
    print("ROBUSTNESS ANALYSIS COMPLETE")
    print("="*60)

    return ts_df, tau_df, fit_df


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Robustness analysis for IMP aging simulation")
    p.add_argument("--csv", type=str, required=True, help="Input CSV from imp_aging.py")
    p.add_argument("--outdir", type=str, default="analysis/robustness", help="Output directory")

    args = p.parse_args()

    run_robustness_analysis(args.csv, args.outdir)
