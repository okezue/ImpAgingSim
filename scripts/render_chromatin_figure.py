#!/usr/bin/env python3
"""Render chromatin-proxy campaign results from committed numerical summaries.

Input: analysis/chromatin_memory/{per_condition,per_run}.csv.
All points/cells are measured condition means; error bars are seed SEM (n=4).
Cluster and mark-fraction means use the trailing 6,250 tau of each 12,500 tau run.
Panel C shows measured first 1/e crossing times, not invented correlation curves.
The site-mark observable is time-centered per site. Its finite-window estimate
should not be equated to the infinite-duration kinetic expectation.
Panel D plots every condition without fitting a universal ratio-only law.

This is a coarse-grained dynamic-label polymer model inspired by chromatin.
No explicit reader/writer proteins, replication, or biological time calibration
is represented in this figure.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analysis/chromatin_memory"
OUT = ROOT / "docs/figures"
INK = "#243345"
MUTED = "#657385"
TEAL = "#007F86"
ORANGE = "#D97943"
PURPLE = "#7556A5"
GRID = "#E4E8ED"


def table(name: str) -> list[dict[str, float | str]]:
    with (SOURCE / name).open(newline="") as handle:
        result = []
        for row in csv.DictReader(handle):
            result.append({
                key: (value if key == "run_id" else float(value) if value else np.nan)
                for key, value in row.items()
            })
    return result


def select(rows, **values):
    return [row for row in rows if all(np.isclose(row[key], value) for key, value in values.items())]


def verify_summaries(conditions, runs):
    """Verify displayed condition summaries against the supplied seed-level data."""
    assert len(conditions) == 40 and len(runs) == 160
    for row in conditions:
        group = select(runs, **{key: row[key] for key in ("eps_BB", "k_off", "k_fb")})
        assert len(group) == int(row["n_seeds"]) == 4
        fields = ["f_B", "largest_B_cluster_fraction"]
        if row["k_fb"] == 0:
            fields += ["tau_BB_peak", "mark_memory_time"]
        for field in fields:
            values = np.array([seed[field] for seed in group])
            assert np.isfinite(values).all()
            np.testing.assert_allclose(values.mean(), row[field + "_mean"], rtol=1e-12)
            np.testing.assert_allclose(values.std(ddof=1) / 2, row[field + "_sem"], rtol=1e-10, atol=1e-13)


def axis_style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for edge in ("left", "bottom"):
        ax.spines[edge].set_color(MUTED)
        ax.spines[edge].set_linewidth(0.8)
    ax.tick_params(length=3.5, width=0.8, color=MUTED, labelsize=10.5, pad=5)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)


def main():
    conditions, runs = table("per_condition.csv"), table("per_run.csv")
    verify_summaries(conditions, runs)
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "svg.fonttype": "path",
        "svg.hashsalt": "ImpAgingSim-chromatin-v1", "mathtext.fontset": "dejavusans",
    })
    fig = plt.figure(figsize=(12, 7.7))
    a = fig.add_axes((0.085, 0.640, 0.335, 0.280))
    b = fig.add_axes((0.565, 0.640, 0.335, 0.280))
    c = fig.add_axes((0.085, 0.125, 0.355, 0.365))
    d = fig.add_axes((0.565, 0.125, 0.355, 0.365))

    for label, x, y, title in (
        ("A", 0.030, 0.953, r"$\epsilon_{BB}=1.0$"),
        ("B", 0.510, 0.953, r"$\epsilon_{BB}=1.5$"),
        ("C", 0.030, 0.515, ""),
        ("D", 0.510, 0.515, ""),
    ):
        fig.text(x, y, label, fontsize=17, weight="bold")
        fig.text(x + 0.055, y, title, fontsize=13)

    koff = sorted({row["k_off"] for row in conditions})
    kfb = sorted({row["k_fb"] for row in conditions})
    cmap = LinearSegmentedColormap.from_list("marked_connectivity", ["#F4F0FA", "#C6B5DD", PURPLE])
    for ax, eps in ((a, 1.0), (b, 1.5)):
        values = np.array([
            [select(conditions, eps_BB=eps, k_off=off, k_fb=fb)[0]["largest_B_cluster_fraction_mean"]
             for fb in kfb] for off in koff
        ])
        heat = ax.pcolormesh(np.arange(len(kfb) + 1) - 0.5, np.arange(len(koff) + 1) - 0.5,
                             values, vmin=0, vmax=1, cmap=cmap, shading="flat",
                             edgecolors="white", linewidth=1.2, rasterized=False)
        ax.set_xticks(range(len(kfb)), ["0", "0.01", "0.03", "0.1", "0.3"])
        ax.set_yticks(range(len(koff)), ["0.001", "0.003", "0.01", "0.03"])
        ax.set_xticks(np.arange(-0.5, len(kfb), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(koff), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.5)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.tick_params(which="major", length=0, pad=6, labelsize=10.5)
        ax.set_xlabel(r"Feedback $k_{\rm fb}$ ($\tau^{-1}$)", labelpad=8)
        ax.set_ylabel(r"Turnover $k_{\rm off}$ ($\tau^{-1}$)", labelpad=8)
        for spine in ax.spines.values():
            spine.set_visible(False)
    colorbar_ax = fig.add_axes((0.922, 0.640, 0.014, 0.280))
    colorbar = fig.colorbar(heat, cax=colorbar_ax, ticks=[0, 0.5, 1])
    colorbar.solids.set_rasterized(False)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=0, labelsize=10)
    colorbar.set_label("Largest-cluster fraction", fontsize=10.5, color=MUTED, labelpad=8)

    axis_style(c)
    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xlim(0.00078, 0.039)
    c.set_ylim(9, 950)
    c.xaxis.set_major_locator(FixedLocator(koff))
    c.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    c.xaxis.set_minor_locator(NullLocator())
    c.yaxis.set_major_locator(FixedLocator([10, 100, 1000]))
    c.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    c.yaxis.set_minor_locator(NullLocator())
    for eps, color, marker in ((1.0, TEAL, "o"), (1.5, ORANGE, "s")):
        rows = sorted(select(conditions, eps_BB=eps, k_fb=0), key=lambda r: r["k_off"])
        c.errorbar([r["k_off"] for r in rows], [r["tau_BB_peak_mean"] for r in rows],
                   yerr=[r["tau_BB_peak_sem"] for r in rows], color=color, marker=marker,
                   markersize=5.5, linewidth=1.8, capsize=3, elinewidth=1.0,
                   markeredgecolor="white", markeredgewidth=0.6,
                   label=rf"B density, $\epsilon_{{BB}}={eps:g}$")
    rows = sorted(select(conditions, eps_BB=1.0, k_fb=0), key=lambda r: r["k_off"])
    c.errorbar([r["k_off"] for r in rows], [r["mark_memory_time_mean"] for r in rows],
               yerr=[r["mark_memory_time_sem"] for r in rows], color=INK, marker="D",
               linestyle=(0, (3, 2)), markersize=4.5, linewidth=1.5, capsize=3,
               markeredgecolor="white", markeredgewidth=0.6, label="Site marks")
    c.set_xlabel(r"Turnover $k_{\rm off}$ ($\tau^{-1}$)", labelpad=8)
    c.set_ylabel(r"Measured $1/e$ time ($\tau$)", labelpad=8)
    c.legend(loc="upper right", frameon=False, fontsize=9.5, handlelength=2.2, labelspacing=0.4)
    c.text(0.04, 0.065, r"$k_{\rm fb}=0$", transform=c.transAxes, color=MUTED, size=10.5)

    axis_style(d)
    d.set_xscale("symlog", linthresh=0.3, linscale=0.55)
    d.set_xlim(-0.025, 420)
    d.set_ylim(0.255, 1.035)
    d.xaxis.set_major_locator(FixedLocator([0, 1, 10, 100, 300]))
    d.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    d.xaxis.set_minor_locator(NullLocator())
    d.set_yticks([0.3, 0.5, 0.7, 0.9, 1.0])
    d.axhline(0.3, color=MUTED, linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
    for eps, color, marker in ((1.0, TEAL, "o"), (1.5, ORANGE, "s")):
        rows = select(conditions, eps_BB=eps)
        d.errorbar([r["k_fb"] / r["k_off"] for r in rows], [r["f_B_mean"] for r in rows],
                   yerr=[r["f_B_sem"] for r in rows], color=color, marker=marker,
                   markersize=5.0, linestyle="none", capsize=2.5, elinewidth=1.0,
                   markeredgecolor="white", markeredgewidth=0.6, alpha=0.9,
                   label=rf"$\epsilon_{{BB}}={eps:g}$")
    d.set_xlabel(r"Feedback / turnover, $k_{\rm fb}/k_{\rm off}$", labelpad=8)
    d.set_ylabel(r"Marked fraction $\langle f_B\rangle$", labelpad=8)
    d.legend(loc="upper left", frameon=False, fontsize=10.5, handlelength=1.5)

    description = (
        "Dynamic-label polymer proxy inspired by chromatin: 160 simulations, four seeds "
        "per condition. A and B: fraction of B beads in the largest contact cluster. "
        "C: measured 1/e times of B-density and site-mark correlations without feedback. "
        "D: trailing-half mean marked fraction versus the feedback/turnover ratio. "
        "Error bars are across-seed SEM; no raw correlation curves or biological claims "
        "are synthesized. Sources: analysis/chromatin_memory/per_condition.csv and per_run.csv."
    )
    fig.savefig(OUT / "chromatin-memory-results.svg", metadata={"Date": None, "Description": description})
    fig.savefig(OUT / "chromatin-memory-results.png", dpi=180, metadata={"Description": description})
    plt.close(fig)
    print(OUT / "chromatin-memory-results.svg")
    print(OUT / "chromatin-memory-results.png")


if __name__ == "__main__":
    main()
