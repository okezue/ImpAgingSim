#!/usr/bin/env python3
"""Render the README's measured static-structure figure from committed data.

Run from any directory with Python, NumPy, and Matplotlib installed:
    python scripts/render_structure_figure.py

Source: the corrected fixed_density_pi099_v2 campaign's exact reciprocal-shell
analysis. Each condition has five independent seeds. Each seed spectrum is the
mean of the final five recorded configurations (steps 242000--250000). Error
bars/bands are the across-seed standard error, as supplied in the source tables.
The displayed channel is S_psi,psi^(N) / 2, where S_psi,psi^(N) uses total-bead
normalization. This is not the density-orthogonal Bhatia--Thornton S_cc.

Only the line segments joining measured shell values are interpolated visually;
there is no smoothing, synthetic morphology, fit, or extrapolation below k_min.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "output/melt/fixed_density_size/fixed_density_pi099_v2/analysis"
OUT = ROOT / "docs/figures"
INK = "#243345"
MUTED = "#657385"
PURPLE = "#7556A5"
GRID = "#E4E8ED"


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def select(table: list[dict[str, str]], **conditions: int) -> list[dict[str, str]]:
    return [r for r in table if all(float(r[key]) == value for key, value in conditions.items())]


def values(table: list[dict[str, str]], field: str) -> np.ndarray:
    return np.array([float(r[field]) for r in table])


def style_axis(ax: plt.Axes, lower: float, upper: float) -> None:
    ax.set_yscale("log")
    ax.set_ylim(lower, upper)
    ax.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    for edge in ("top", "right"):
        ax.spines[edge].set_visible(False)
    for edge in ("bottom", "left"):
        ax.spines[edge].set_color(MUTED)
        ax.spines[edge].set_linewidth(0.8)
    ax.tick_params(length=3.5, width=0.8, color=MUTED, pad=6)


def main() -> None:
    spectra = rows("shell_spectra_condition.csv")
    summary = rows("condition_summary.csv")
    assert len(summary) == 6 and {int(r["n_seeds"]) for r in summary} == {5}
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.labelsize": 13,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "svg.fonttype": "path",
            "svg.hashsalt": "ImpAgingSim-measured-structure-v1",
            "mathtext.fontset": "dejavusans",
        }
    )
    fig = plt.figure(figsize=(12, 5.2), dpi=150)
    ax_a = fig.add_axes([0.09, 0.19, 0.425, 0.63])
    ax_b = fig.add_axes([0.67, 0.19, 0.29, 0.63])

    fig.text(0.025, 0.90, "A", weight="bold", fontsize=17)
    fig.text(0.09, 0.90, "Measured composition spectrum", weight="bold", fontsize=14)
    fig.text(0.09, 0.852, r"$M=144$ chains", color=MUTED, fontsize=11)
    fig.text(0.60, 0.90, "B", weight="bold", fontsize=17)
    fig.text(0.67, 0.90, "System-size comparison", weight="bold", fontsize=14)
    fig.text(0.67, 0.852, r"Fixed bead density", color=MUTED, fontsize=11)

    style_axis(ax_a, 0.18, 1100)
    style_axis(ax_b, 3, 1300)
    ax_a.set_xlim(0.20, 1.53)
    ax_a.set_xticks([0.3, 0.6, 0.9, 1.2, 1.5])
    ax_a.set_xlabel(r"Wavevector, $k\sigma$", labelpad=9)
    ax_a.set_ylabel(r"$S_{\psi\psi}^{(N)}(k)\,/\,2$", labelpad=12)
    ax_b.set_xlim(100, 620)
    ax_b.set_xticks([144, 288, 576])
    ax_b.set_xlabel(r"Number of chains, $M$", labelpad=9)
    ax_b.set_ylabel(r"Peak $S_{\psi\psi}^{(N)}\,/\,2$", labelpad=10)

    for kappa, color in ((0, MUTED), (1, PURPLE)):
        s = sorted(select(spectra, n_chains=144, kappa=kappa), key=lambda r: float(r["q"]))
        assert len(s) == 24
        q = values(s, "q")
        mean = values(s, "S_psi_psi_over_2_mean")
        sem = values(s, "S_psi_psi_over_2_sem")
        assert np.all(mean > sem)
        ax_a.fill_between(q, mean - sem, mean + sem, color=color, alpha=0.15, linewidth=0)
        ax_a.plot(q, mean, color=color, linewidth=2.0, marker="o", markersize=3.5,
                  markeredgewidth=0, label=rf"$\kappa={kappa}$")

        c = sorted(select(summary, kappa=kappa), key=lambda r: int(r["n_chains"]))
        x = values(c, "n_chains")
        y = values(c, "primary_amplitude_seed_mean_at_selected_shell")
        error = values(c, "primary_amplitude_seed_sem_at_selected_shell")
        ax_b.errorbar(x, y, yerr=error, color=color, linewidth=2, marker="o",
                      markersize=6.5, markeredgecolor="white", markeredgewidth=0.7,
                      capsize=4, elinewidth=1.3, capthick=1.3)

    reference = select(summary, n_chains=144, kappa=1)[0]
    kmin = 2 * np.pi / float(reference["box_size"])
    maximum = float(reference["primary_amplitude_seed_mean_at_selected_shell"])
    ax_a.axvline(kmin, color=MUTED, linewidth=0.9, linestyle=(0, (3, 3)), zorder=0)
    ax_a.annotate("Box-limited maximum", xy=(kmin, maximum),
                  xytext=(0.62, 460), ha="left", va="center", fontsize=10.5,
                  color=PURPLE, arrowprops={"arrowstyle": "-", "color": PURPLE,
                                            "linewidth": 0.9, "shrinkB": 7})
    ax_a.text(kmin + 0.022, 0.24, r"$k_{\min}=2\pi/L$", fontsize=10, color=MUTED)
    # Direct labels preserve the two-channel key without occupying a data region.
    last_correlated = select(summary, n_chains=576, kappa=1)[0]
    last_correlated_peak = float(last_correlated["primary_amplitude_seed_mean_at_selected_shell"])
    ax_b.text(606, last_correlated_peak, r"$\kappa=1$", color=PURPLE,
              fontsize=12, ha="right", va="bottom")
    ax_b.text(606, 8.45, r"$\kappa=0$", color=MUTED,
              fontsize=12, ha="right", va="bottom")
    fig.text(0.09, 0.045, r"$N_{\rm chain}=40$   ·   $f_A=0.5$   ·   $\pi=0.99$   ·   $T^*=0.7$",
             fontsize=11, color=MUTED)
    fig.text(0.96, 0.045, "Mean ± SEM · 5 independent seeds", ha="right", fontsize=11, color=MUTED)

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        options = {"metadata": {"Date": None}} if ext == "svg" else {}
        fig.savefig(OUT / f"measured-structure.{ext}", dpi=180, **options)
    plt.close(fig)


if __name__ == "__main__":
    main()
