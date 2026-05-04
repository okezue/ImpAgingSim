#!/usr/bin/env python3
"""Create comparison overlay figures for iid vs correlated ensembles."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
matplotlib.rcParams.update({
    'font.size':11,'axes.labelsize':13,'axes.titlesize':13,
    'legend.fontsize':9,'figure.dpi':200,'savefig.dpi':200,
    'lines.linewidth':1.5,'lines.markersize':5,
    'axes.grid':True,'grid.alpha':0.3,
})
def main():
    df = pd.read_csv("output/results.csv")

    outdir = "analysis/comparison"
    os.makedirs(outdir, exist_ok=True)

    epsilons = sorted(df["epsilon"].unique())
    tws = sorted(df["tw"].unique())

    for eps in epsilons:
        # Figure 1: Q_mean comparison (normalized by Q at lag=min)
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Left panel: Q_norm (collapsed) for both ensembles
        ax = axes[0]
        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            sub = df[(df["ensemble"] == ens) & (df["epsilon"] == eps)].copy()
            gcols = ["tw", "lag"]
            agg = sub.groupby(gcols).agg(Q_mean=("Q_mean", "mean")).reset_index()

            for tw in tws:
                tw_data = agg[agg["tw"] == tw].sort_values("lag")
                if len(tw_data) > 0:
                    Q0 = tw_data["Q_mean"].iloc[0]
                    if Q0 > 0:
                        Q_norm = tw_data["Q_mean"] / Q0
                        alpha = 0.3 + 0.7 * (tws.index(tw) / max(1, len(tws)-1))
                        ax.plot(tw_data["lag"], Q_norm, marker=marker, linestyle="-",
                               color=color, alpha=alpha, markersize=4,
                               label=f"{ens} tw={int(tw)}" if tw == tws[0] else None)

        ax.axhline(y=np.exp(-1), color='gray', linestyle='--', alpha=0.5, label=r'$e^{-1}$')
        ax.set_xscale("log")
        ax.set_xlabel("lag (MC sweeps)")
        ax.set_ylabel(r"$Q_{norm} = Q(t_w, t) / Q(t_w, 0)$")
        ax.set_title(f"Q collapse comparison | epsilon={eps}")
        ax.legend(loc="best", fontsize=8)
        ax.set_ylim(0, 1.1)

        # Right panel: chi4 for both ensembles
        ax = axes[1]
        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            sub = df[(df["ensemble"] == ens) & (df["epsilon"] == eps)].copy()
            gcols = ["tw", "lag"]
            agg = sub.groupby(gcols).agg(chi4=("chi4", "mean")).reset_index()

            for tw in tws:
                tw_data = agg[agg["tw"] == tw].sort_values("lag")
                if len(tw_data) > 0:
                    alpha = 0.3 + 0.7 * (tws.index(tw) / max(1, len(tws)-1))
                    ax.plot(tw_data["lag"], tw_data["chi4"], marker=marker, linestyle="-",
                           color=color, alpha=alpha, markersize=4,
                           label=f"{ens} tw={int(tw)}" if tw == tws[0] else None)

        ax.set_xscale("log")
        ax.set_xlabel("lag (MC sweeps)")
        ax.set_ylabel(r"$\chi_4$")
        ax.set_title(f"chi4 comparison | epsilon={eps}")
        ax.legend(loc="best", fontsize=8)

        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"Fig_compare_collapses_eps-{int(eps)}.png"), dpi=200)
        plt.savefig(os.path.join(outdir, f"Fig_compare_collapses_eps-{int(eps)}.pdf"))
        plt.close()
        print(f"Saved comparison figure for epsilon={eps}")

    # Create summary figure: chi4_star vs tw for all conditions
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, eps in enumerate(epsilons):
        ax = axes[i]
        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            sub = df[(df["ensemble"] == ens) & (df["epsilon"] == eps)].copy()
            gcols = ["tw", "lag"]
            agg = sub.groupby(gcols).agg(chi4=("chi4", "mean")).reset_index()

            chi4_stars = []
            for tw in tws:
                tw_data = agg[agg["tw"] == tw]
                if len(tw_data) > 0:
                    chi4_stars.append((tw, tw_data["chi4"].max()))

            if chi4_stars:
                tws_plot, chi4_plot = zip(*chi4_stars)
                ax.plot(tws_plot, chi4_plot, marker=marker, linestyle="-",
                       color=color, label=ens)

        ax.set_xscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)")
        ax.set_ylabel(r"$\chi_4^*$")
        ax.set_title(f"epsilon={int(eps)}")
        ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_chi4_star_vs_tw.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_chi4_star_vs_tw.pdf"))
    plt.close()
    print("Saved chi4_star vs tw summary figure")

    # Create t_star vs tw figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, eps in enumerate(epsilons):
        ax = axes[i]
        for ens, color, marker in [("iid", "blue", "o"), ("correlated", "red", "s")]:
            sub = df[(df["ensemble"] == ens) & (df["epsilon"] == eps)].copy()
            gcols = ["tw", "lag"]
            agg = sub.groupby(gcols).agg(chi4=("chi4", "mean")).reset_index()

            t_stars = []
            for tw in tws:
                tw_data = agg[agg["tw"] == tw].sort_values("lag")
                if len(tw_data) > 0:
                    idx_max = tw_data["chi4"].idxmax()
                    t_star = tw_data.loc[idx_max, "lag"]
                    t_stars.append((tw, t_star))

            if t_stars:
                tws_plot, tstar_plot = zip(*t_stars)
                ax.plot(tws_plot, tstar_plot, marker=marker, linestyle="-",
                       color=color, label=ens)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$t_w$ (MC sweeps)")
        ax.set_ylabel(r"$t^*$ (lag at $\chi_4$ peak)")
        ax.set_title(f"epsilon={int(eps)}")
        ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "Fig_tstar_vs_tw.png"), dpi=200)
    plt.savefig(os.path.join(outdir, "Fig_tstar_vs_tw.pdf"))
    plt.close()
    print("Saved t_star vs tw summary figure")

if __name__ == "__main__":
    main()
