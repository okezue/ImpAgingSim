#!/usr/bin/env python3
"""Reproduce the README protein-condensation figure from committed source tables.

Run from any directory with NumPy and Matplotlib installed:
    python scripts/render_protein_figure.py

Source: analysis/protein_condensation/{per_run,per_condition}.csv. All 360 runs
contribute to the two heatmaps (90 conditions, four seeds each). The lower panels
show the f_B=0.3 slice at kappa=0, 0.5, and 1; all nine sampled attractions are
included. Source means and SEMs are verified against per-run summaries before
plotting. The raw trajectories are not required or represented by this figure.

Each per-run observable is averaged over steps >500000 through 1000000, i.e.
250 saved configurations at 2000-step intervals, with dt=0.005 tau. Error bars
are SEM across four seeds, not temporal fluctuations or reciprocal-shell error.
B contacts use a periodic distance <=1.5 sigma, including bonded neighbors.
The contact graph alone does not identify individual chains inside a cluster.
No equilibrium phase boundary, morphology, kinetic mechanism, or unobserved
parameter interpolation is inferred. Heatmap cells represent sampled conditions;
their edges are the midpoints between sampled attraction values.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "analysis/protein_condensation"
OUT = ROOT / "docs/figures"
INK = "#243345"
MUTED = "#657385"
PURPLE = "#7556A5"
MIDPURPLE = "#B0A0C8"
RULE = "#E4E8ED"
OBSERVABLES = ("largest_B_cluster_fraction", "mean_Rg", "mean_BB_coordination")
KEYS = ("eps_BB", "kappa", "f_A")


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def pick(table: list[dict[str, str]], **filters: float) -> list[dict[str, str]]:
    return [r for r in table if all(float(r[k]) == v for k, v in filters.items())]


def array(rows: list[dict[str, str]], key: str) -> np.ndarray:
    return np.array([float(r[key]) for r in rows])


def verified_conditions() -> list[dict[str, str]]:
    conditions = read_csv("per_condition.csv")
    runs = read_csv("per_run.csv")
    assert len(conditions) == 90 and len(runs) == 360
    for condition in conditions:
        group = pick(runs, **{key: float(condition[key]) for key in KEYS})
        assert len(group) == int(condition["n_seeds"]) == 4
        for observable in OBSERVABLES:
            values = array(group, observable)
            assert np.all(np.isfinite(values))
            np.testing.assert_allclose(values.mean(), float(condition[f"{observable}_mean"]), rtol=1e-12)
            np.testing.assert_allclose(values.std(ddof=1) / np.sqrt(4),
                                       float(condition[f"{observable}_sem"]), rtol=1e-12)
    return conditions


def bin_edges(centers: np.ndarray) -> np.ndarray:
    return np.r_[centers[0] - (centers[1] - centers[0]) / 2,
                 (centers[:-1] + centers[1:]) / 2,
                 centers[-1] + (centers[-1] - centers[-2]) / 2]


def style_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    for edge in ("bottom", "left"):
        ax.spines[edge].set_color(MUTED)
        ax.spines[edge].set_linewidth(0.8)
    ax.tick_params(length=3.5, width=0.8, color=MUTED, pad=6)
    ax.set_axisbelow(True)


def main() -> None:
    table = verified_conditions()
    eps = np.array(sorted({float(c["eps_BB"]) for c in table}))
    kappas = np.array(sorted({float(c["kappa"]) for c in table}))
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 12,
        "axes.labelsize": 13, "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "svg.fonttype": "path",
        "svg.hashsalt": "ImpAgingSim-protein-condensation-v1",
        "mathtext.fontset": "dejavusans",
    })
    fig = plt.figure(figsize=(12, 7.8), dpi=150)
    heat_axes = [fig.add_axes([x, 0.60, 0.345, 0.30]) for x in (0.095, 0.555)]
    line_axes = [fig.add_axes([x, 0.11, 0.345, 0.32]) for x in (0.095, 0.555)]
    cmap = LinearSegmentedColormap.from_list("cluster", ["#F4F2F8", "#DED7EC", "#A18ABF", PURPLE, "#36264F"])
    norm = Normalize(vmin=0, vmax=1)

    for ax, f_A, letter in zip(heat_axes, (0.7, 0.5), ("A", "B")):
        grid = np.empty((len(kappas), len(eps)))
        for i, kappa in enumerate(kappas):
            for j, attraction in enumerate(eps):
                row = pick(table, f_A=f_A, kappa=kappa, eps_BB=attraction)
                assert len(row) == 1
                grid[i, j] = float(row[0]["largest_B_cluster_fraction_mean"])
        image = ax.pcolormesh(bin_edges(eps), bin_edges(kappas), grid, cmap=cmap,
                             norm=norm, edgecolors="white", linewidth=1.2, shading="flat")
        ax.set_xticks(eps, [f"{x:g}" for x in eps])
        ax.tick_params(axis="x", labelsize=10.5)
        ax.set_yticks(kappas, [f"{x:g}" for x in kappas])
        ax.set_xlabel(r"B–B attraction, $\epsilon_{BB}/\epsilon_{\rm core}$", labelpad=9)
        ax.set_ylabel(r"Sequence correlation, $\kappa$", labelpad=9)
        style_axis(ax)
        x = ax.get_position().x0
        fig.text(x - 0.055, 0.926, letter, fontsize=17, weight="bold")
        fig.text(x, 0.926, rf"$f_B={1-f_A:.1f}$", fontsize=13, color=MUTED)
    cax = fig.add_axes([0.929, 0.60, 0.014, 0.30])
    colorbar = fig.colorbar(image, cax=cax, ticks=[0, 0.5, 1])
    colorbar.solids.set_rasterized(False)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=2.5, labelsize=11, pad=4)
    colorbar.set_label("Largest-cluster fraction", rotation=90, fontsize=11.5, labelpad=5)

    styles = ((0.0, MUTED, "-", "o"), (0.5, MIDPURPLE, "--", "s"), (1.0, PURPLE, "-", "o"))
    for ax, metric, letter, ylabel in zip(
        line_axes, ("mean_Rg", "mean_BB_coordination"), ("C", "D"),
        (r"Mean $R_g/\sigma$", "B neighbors per B bead"),
    ):
        for kappa, color, linestyle, marker in styles:
            rows = sorted(pick(table, f_A=0.7, kappa=kappa), key=lambda r: float(r["eps_BB"]))
            ax.errorbar(array(rows, "eps_BB"), array(rows, f"{metric}_mean"),
                        yerr=array(rows, f"{metric}_sem"), color=color, linewidth=2.2,
                        linestyle=linestyle, marker=marker, markersize=5.2,
                        markeredgecolor="white", markeredgewidth=0.65,
                        capsize=3, elinewidth=1.1, capthick=1.1, label=rf"$\kappa={kappa:g}$")
        style_axis(ax)
        ax.grid(axis="y", color=RULE, linewidth=0.7)
        ax.set_xlim(0.38, 3.12)
        ax.set_xticks([0.5, 1, 1.5, 2, 2.5, 3])
        ax.set_xlabel(r"B–B attraction, $\epsilon_{BB}/\epsilon_{\rm core}$", labelpad=9)
        ax.set_ylabel(ylabel, labelpad=9)
        x = ax.get_position().x0
        fig.text(x - 0.055, 0.456, letter, fontsize=17, weight="bold")
        fig.text(x, 0.456, r"$f_B=0.3$", fontsize=13, color=MUTED)
    line_axes[0].set_ylim(3.32, 3.91)
    line_axes[0].set_yticks([3.4, 3.6, 3.8])
    line_axes[1].set_ylim(0, 14.8)
    line_axes[1].set_yticks([0, 4, 8, 12])
    line_axes[1].legend(loc="upper left", frameon=False, fontsize=11, ncol=3,
                        handlelength=1.4, handletextpad=0.4, columnspacing=0.8, borderaxespad=0.2)
    OUT.mkdir(parents=True, exist_ok=True)
    description = (
        "Measured largest B-cluster fraction in the 360-run protein-condensation campaign, "
        "with chain size and B-B coordination at f_B=0.3. Conditions are averaged over "
        "four seeds. Per-run values use the trailing half of 1,000,000 production steps. "
        "No equilibrium phase boundary is inferred."
    )
    fig.savefig(OUT / "protein-condensation-results.svg", metadata={"Date": None, "Description": description})
    fig.savefig(OUT / "protein-condensation-results.png", dpi=180, metadata={"Description": description})
    plt.close(fig)
    print(f"Verified 90 conditions against 360 per-run summaries; wrote {OUT / 'protein-condensation-results.svg'}")


if __name__ == "__main__":
    main()
