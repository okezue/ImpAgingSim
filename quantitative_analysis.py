#!/usr/bin/env python3
"""Quantitative analysis: power-law fits and ensemble separation metrics."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress

def power_law(x, a, b):
    """Power law: y = a * x^b"""
    return a * np.power(x, b)

def log_power_law(log_x, log_a, b):
    """Linear in log-log space: log(y) = log(a) + b*log(x)"""
    return log_a + b * log_x

def extract_timescales(df):
    """Extract t* and chi4* for each (ensemble, epsilon, tw)."""
    results = []

    for (ens, eps), group in df.groupby(["ensemble", "epsilon"]):
        # Aggregate over disorder samples
        agg = group.groupby(["tw", "lag"]).agg(
            chi4_mean=("chi4", "mean"),
            chi4_std=("chi4", "std"),
            Q_mean=("Q_mean", "mean"),
        ).reset_index()

        for tw in agg["tw"].unique():
            tw_data = agg[agg["tw"] == tw].sort_values("lag")
            if len(tw_data) < 3:
                continue

            # Find t* (lag at chi4 peak)
            idx_max = tw_data["chi4_mean"].idxmax()
            t_star = tw_data.loc[idx_max, "lag"]
            chi4_star = tw_data.loc[idx_max, "chi4_mean"]
            chi4_star_std = tw_data.loc[idx_max, "chi4_std"]

            # Compute Q_norm and find tau_alpha if possible
            Q0 = tw_data["Q_mean"].iloc[0]
            if Q0 > 0:
                tw_data = tw_data.copy()
                tw_data["Q_norm"] = tw_data["Q_mean"] / Q0
                # Find crossing through e^-1
                threshold = np.exp(-1)
                below = tw_data[tw_data["Q_norm"] < threshold]
                if len(below) > 0:
                    tau_alpha = below["lag"].iloc[0]
                else:
                    tau_alpha = np.nan
            else:
                tau_alpha = np.nan

            results.append({
                "ensemble": ens,
                "epsilon": eps,
                "tw": tw,
                "t_star": t_star,
                "chi4_star": chi4_star,
                "chi4_star_std": chi4_star_std,
                "tau_alpha": tau_alpha,
            })

    return pd.DataFrame(results)

def fit_power_law(x, y, min_points=3):
    """Fit power law y = a * x^b using log-log linear regression."""
    # Filter out zeros and invalid values
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    x_valid = x[mask]
    y_valid = y[mask]

    if len(x_valid) < min_points:
        return None, None, None, None

    log_x = np.log10(x_valid)
    log_y = np.log10(y_valid)

    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)

    a = 10**intercept
    b = slope
    r_squared = r_value**2

    return a, b, r_squared, std_err

def main():
    df = pd.read_csv("output/results.csv")
    outdir = "analysis/quantitative"
    os.makedirs(outdir, exist_ok=True)

    # Extract timescales
    timescales = extract_timescales(df)
    timescales.to_csv(os.path.join(outdir, "timescales.csv"), index=False)
    print(f"Saved timescales to {outdir}/timescales.csv")
    print(timescales.to_string())
    print()

    # Filter to tw > 0 for power law fitting
    ts_nonzero = timescales[timescales["tw"] > 0].copy()

    epsilons = sorted(timescales["epsilon"].unique())
    ensembles = ["iid", "correlated"]

    # =========================================================================
    # 1. Power-law fits for t*(tw)
    # =========================================================================
    print("=" * 60)
    print("POWER-LAW FITS: t*(t_w) ~ A * t_w^nu")
    print("=" * 60)

    fit_results = []

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for i, eps in enumerate(epsilons):
        ax = axes[i]

        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            subset = ts_nonzero[(ts_nonzero["ensemble"] == ens) &
                                (ts_nonzero["epsilon"] == eps)]

            if len(subset) < 2:
                continue

            tw = subset["tw"].values
            t_star = subset["t_star"].values

            # Fit power law
            a, b, r2, se = fit_power_law(tw, t_star)

            if a is not None:
                fit_results.append({
                    "quantity": "t_star",
                    "ensemble": ens,
                    "epsilon": eps,
                    "A": a,
                    "exponent": b,
                    "R_squared": r2,
                    "std_err": se,
                })

                print(f"  eps={eps}, {ens}: t* ~ {a:.4f} * t_w^{b:.4f}  (R²={r2:.4f})")

                # Plot fit line
                tw_fit = np.logspace(np.log10(tw.min()), np.log10(tw.max()), 50)
                t_star_fit = power_law(tw_fit, a, b)
                ax.plot(tw_fit, t_star_fit, '--', color=color, alpha=0.7,
                       label=f"{ens}: $\\nu$={b:.3f}")

            # Plot data
            ax.scatter(tw, t_star, color=color, marker=marker, s=60, zorder=5)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=11)
        ax.set_ylabel(r"$t^*$ (lag at $\chi_4$ peak)", fontsize=11)
        ax.set_title(f"$\\epsilon$ = {int(eps)}", fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_tstar_powerlaw_fits.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_tstar_powerlaw_fits.pdf"))
    plt.close()
    print(f"\nSaved: {outdir}/Fig_tstar_powerlaw_fits.png")

    # =========================================================================
    # 2. Power-law fits for chi4*(tw)
    # =========================================================================
    print("\n" + "=" * 60)
    print("POWER-LAW FITS: chi4*(t_w) ~ A * t_w^mu (early-time growth)")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for i, eps in enumerate(epsilons):
        ax = axes[i]

        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            subset = ts_nonzero[(ts_nonzero["ensemble"] == ens) &
                                (ts_nonzero["epsilon"] == eps)]

            if len(subset) < 2:
                continue

            tw = subset["tw"].values
            chi4_star = subset["chi4_star"].values

            # Fit only early-time growth (before peak) - use first 3-4 points
            n_fit = min(4, len(tw))
            a, b, r2, se = fit_power_law(tw[:n_fit], chi4_star[:n_fit])

            if a is not None:
                fit_results.append({
                    "quantity": "chi4_star_early",
                    "ensemble": ens,
                    "epsilon": eps,
                    "A": a,
                    "exponent": b,
                    "R_squared": r2,
                    "std_err": se,
                })

                print(f"  eps={eps}, {ens}: chi4* ~ {a:.4f} * t_w^{b:.4f}  (R²={r2:.4f}, n={n_fit})")

            # Plot data with error bars
            ax.errorbar(tw, chi4_star, yerr=subset["chi4_star_std"].values,
                       color=color, marker=marker, linestyle='-', capsize=3,
                       label=f"{ens}")

        ax.set_xscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=11)
        ax.set_ylabel(r"$\chi_4^*$", fontsize=11)
        ax.set_title(f"$\\epsilon$ = {int(eps)}", fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_chi4star_vs_tw.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_chi4star_vs_tw.pdf"))
    plt.close()
    print(f"\nSaved: {outdir}/Fig_chi4star_vs_tw.png")

    # =========================================================================
    # 3. Ensemble separation metrics
    # =========================================================================
    print("\n" + "=" * 60)
    print("ENSEMBLE SEPARATION: (iid - correlated) / correlated")
    print("=" * 60)

    separation_results = []

    for eps in epsilons:
        if eps == 0:
            continue  # Skip eps=0 where they're identical

        iid_data = ts_nonzero[(ts_nonzero["ensemble"] == "iid") &
                              (ts_nonzero["epsilon"] == eps)]
        corr_data = ts_nonzero[(ts_nonzero["ensemble"] == "correlated") &
                               (ts_nonzero["epsilon"] == eps)]

        # Merge on tw
        merged = pd.merge(iid_data, corr_data, on="tw", suffixes=("_iid", "_corr"))

        for _, row in merged.iterrows():
            tw = row["tw"]

            # t* separation
            t_star_iid = row["t_star_iid"]
            t_star_corr = row["t_star_corr"]
            if t_star_corr > 0:
                delta_t_star = (t_star_iid - t_star_corr) / t_star_corr
            else:
                delta_t_star = np.nan

            # chi4* separation
            chi4_iid = row["chi4_star_iid"]
            chi4_corr = row["chi4_star_corr"]
            if chi4_corr > 0:
                delta_chi4 = (chi4_iid - chi4_corr) / chi4_corr
            else:
                delta_chi4 = np.nan

            separation_results.append({
                "epsilon": eps,
                "tw": tw,
                "t_star_iid": t_star_iid,
                "t_star_corr": t_star_corr,
                "delta_t_star_rel": delta_t_star,
                "chi4_star_iid": chi4_iid,
                "chi4_star_corr": chi4_corr,
                "delta_chi4_star_rel": delta_chi4,
            })

    sep_df = pd.DataFrame(separation_results)
    sep_df.to_csv(os.path.join(outdir, "ensemble_separation.csv"), index=False)
    print(f"\nSaved: {outdir}/ensemble_separation.csv")
    print(sep_df.to_string())

    # Plot ensemble separation
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: delta_t_star vs tw
    ax = axes[0]
    for eps in [3.0, 6.0]:
        subset = sep_df[sep_df["epsilon"] == eps]
        ax.plot(subset["tw"], subset["delta_t_star_rel"] * 100,
               marker='o', linestyle='-', label=f"$\\epsilon$={int(eps)}")
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=11)
    ax.set_ylabel(r"$\Delta t^* / t^*_{corr}$ (%)", fontsize=11)
    ax.set_title("Relative separation in $t^*$\n(positive = iid > correlated)", fontsize=11)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: delta_chi4 vs tw
    ax = axes[1]
    for eps in [3.0, 6.0]:
        subset = sep_df[sep_df["epsilon"] == eps]
        ax.plot(subset["tw"], subset["delta_chi4_star_rel"] * 100,
               marker='s', linestyle='-', label=f"$\\epsilon$={int(eps)}")
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$t_w$ (MC sweeps)", fontsize=11)
    ax.set_ylabel(r"$\Delta \chi_4^* / \chi_{4,corr}^*$ (%)", fontsize=11)
    ax.set_title("Relative separation in $\\chi_4^*$\n(positive = iid > correlated)", fontsize=11)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_ensemble_separation.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_ensemble_separation.pdf"))
    plt.close()
    print(f"\nSaved: {outdir}/Fig_ensemble_separation.png")

    # =========================================================================
    # 4. Summary statistics
    # =========================================================================
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)

    # Average separation across tw
    for eps in [3.0, 6.0]:
        subset = sep_df[sep_df["epsilon"] == eps]
        mean_delta_t = subset["delta_t_star_rel"].mean() * 100
        std_delta_t = subset["delta_t_star_rel"].std() * 100
        mean_delta_chi4 = subset["delta_chi4_star_rel"].mean() * 100
        std_delta_chi4 = subset["delta_chi4_star_rel"].std() * 100

        print(f"\neps={eps}:")
        print(f"  Mean delta_t*:    {mean_delta_t:+.1f}% +/- {std_delta_t:.1f}%")
        print(f"  Mean delta_chi4*: {mean_delta_chi4:+.1f}% +/- {std_delta_chi4:.1f}%")

    # Save fit results
    fit_df = pd.DataFrame(fit_results)
    fit_df.to_csv(os.path.join(outdir, "powerlaw_fits.csv"), index=False)
    print(f"\nSaved: {outdir}/powerlaw_fits.csv")

    # =========================================================================
    # 5. Exponent comparison figure
    # =========================================================================
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Extract t_star exponents
    t_star_fits = fit_df[fit_df["quantity"] == "t_star"]

    # Left: nu (t* exponent) vs epsilon
    ax = axes[0]
    for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
        subset = t_star_fits[t_star_fits["ensemble"] == ens]
        ax.errorbar(subset["epsilon"], subset["exponent"],
                   yerr=subset["std_err"], color=color, marker=marker,
                   linestyle='-', capsize=4, markersize=8, label=ens)
    ax.set_xlabel(r"$\epsilon$ (disorder strength)", fontsize=11)
    ax.set_ylabel(r"$\nu$ (exponent in $t^* \sim t_w^\nu$)", fontsize=11)
    ax.set_title("Scaling exponent for $t^*(t_w)$", fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Right: Exponent difference
    ax = axes[1]
    iid_fits = t_star_fits[t_star_fits["ensemble"] == "iid"].set_index("epsilon")
    corr_fits = t_star_fits[t_star_fits["ensemble"] == "correlated"].set_index("epsilon")

    common_eps = [e for e in epsilons if e in iid_fits.index and e in corr_fits.index]
    delta_nu = []
    delta_nu_err = []
    for eps in common_eps:
        nu_iid = iid_fits.loc[eps, "exponent"]
        nu_corr = corr_fits.loc[eps, "exponent"]
        se_iid = iid_fits.loc[eps, "std_err"]
        se_corr = corr_fits.loc[eps, "std_err"]
        delta_nu.append(nu_iid - nu_corr)
        delta_nu_err.append(np.sqrt(se_iid**2 + se_corr**2))

    ax.errorbar(common_eps, delta_nu, yerr=delta_nu_err,
               color='purple', marker='D', linestyle='-', capsize=4, markersize=8)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel(r"$\epsilon$ (disorder strength)", fontsize=11)
    ax.set_ylabel(r"$\Delta\nu = \nu_{iid} - \nu_{corr}$", fontsize=11)
    ax.set_title("Exponent difference between ensembles", fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_exponent_comparison.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_exponent_comparison.pdf"))
    plt.close()
    print(f"Saved: {outdir}/Fig_exponent_comparison.png")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

if __name__ == "__main__":
    main()
