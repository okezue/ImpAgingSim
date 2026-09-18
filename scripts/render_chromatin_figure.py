#!/usr/bin/env python3
"""Render chromatin-proxy campaign results from committed numerical summaries.

Input: analysis/chromatin_memory/{per_condition,per_run}.csv.
All points/cells are measured condition means; error bars are seed SEM (n=4).
Cluster and mark-fraction means use the trailing 6,250 tau of each 12,500 tau run.
Panel A pairs both attractions within each cell. Panel B separates density
relaxation from site-mark fluctuations and shows unobserved density crossings
as lower bounds. Panel C shows no-feedback first 1/e crossing times.
The site-mark observable is time-centered per site. Its finite-window estimate
should not be equated to the infinite-duration kinetic expectation.
Panel D plots every condition without fitting a universal ratio-only law.

This is a coarse-grained dynamic-label polymer model inspired by chromatin.
No explicit reader/writer proteins, replication, or biological time calibration
is represented in this figure.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.collections import PolyCollection
from matplotlib.cm import ScalarMappable
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
        if row["k_fb"] == 0 or (row["k_off"] == 0.01 and row["eps_BB"] == 1.5):
            fields += ["mark_memory_time"]
            if np.isfinite(row["tau_BB_peak_mean"]):
                fields += ["tau_BB_peak"]
            else:
                assert all(not np.isfinite(seed["tau_BB_peak"]) for seed in group)
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
    a = fig.add_axes((0.085, 0.645, 0.300, 0.280))
    b = fig.add_axes((0.565, 0.645, 0.355, 0.280))
    c = fig.add_axes((0.085, 0.125, 0.355, 0.365))
    d = fig.add_axes((0.565, 0.125, 0.355, 0.365))
    for label, x, y in (("A", 0.030, 0.958), ("B", 0.510, 0.958),
                        ("C", 0.030, 0.515), ("D", 0.510, 0.515)):
        fig.text(x, y, label, fontsize=17, weight="bold")

    koff = sorted({row["k_off"] for row in conditions})
    kfb = sorted({row["k_fb"] for row in conditions})
    cmap = LinearSegmentedColormap.from_list("marked_connectivity", ["#F4F0FA", "#C6B5DD", PURPLE])
    norm = Normalize(0, 1)
    # Paired measured values, not their mean: upper-left epsilon=1,
    # lower-right epsilon=1.5. Each cell is a sampled categorical condition.
    triangles, colors = [], []
    for j, off in enumerate(koff):
        for i, fb in enumerate(kfb):
            lo, hi, bot, top = i - 0.5, i + 0.5, j - 0.5, j + 0.5
            triangles.extend([[(lo, bot), (lo, top), (hi, top)],
                              [(lo, bot), (hi, top), (hi, bot)]])
            for eps in (1.0, 1.5):
                value = select(conditions, eps_BB=eps, k_off=off, k_fb=fb)[0]["largest_B_cluster_fraction_mean"]
                colors.append(cmap(norm(value)))
    a.add_collection(PolyCollection(triangles, facecolors=colors, edgecolors="white", linewidths=0.8))
    a.set_xlim(-0.5, len(kfb) - 0.5)
    a.set_ylim(-0.5, len(koff) - 0.5)
    a.set_xticks(range(len(kfb)), ["0", "0.01", "0.03", "0.1", "0.3"])
    a.set_yticks(range(len(koff)), ["0.001", "0.003", "0.01", "0.03"])
    a.tick_params(length=0, pad=6, labelsize=10.5)
    a.set_xlabel(r"Feedback $k_{\rm fb}$ ($\tau^{-1}$)", labelpad=8)
    a.set_ylabel(r"Turnover $k_{\rm off}$ ($\tau^{-1}$)", labelpad=8)
    for spine in a.spines.values():
        spine.set_visible(False)
    colorbar_ax = fig.add_axes((0.400, 0.645, 0.013, 0.280))
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=colorbar_ax, ticks=[0, 0.5, 1])
    colorbar.solids.set_rasterized(False)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=0, labelsize=10)
    colorbar.set_label("Largest-cluster fraction", fontsize=10, color=MUTED, labelpad=7)

    # Show actual density relaxation times against site-mark times. The unobserved
    # density crossings are censored at the evaluated lag limit, not finite means.
    settings = json.loads((SOURCE / "manifest.json").read_text())["simulation_settings"]
    saved_steps = np.arange(settings["mode_interval"], settings["n_steps"] + 1, settings["mode_interval"])
    n_frames = np.count_nonzero(saved_steps > 0.5 * settings["n_steps"])
    lag_limit = int(0.5 * (n_frames - 1)) * settings["mode_interval"] * settings["dt"]
    assert lag_limit == 3122.5
    rows = sorted(select(conditions, eps_BB=1.5, k_off=0.01), key=lambda row: row["k_fb"])
    axis_style(b)
    b.set_xscale("symlog", linthresh=0.005, linscale=0.6)
    b.set_yscale("log")
    b.set_xlim(-0.0005, 0.44)
    b.set_ylim(4, 6500)
    b.set_xticks(kfb, ["0", "0.01", "0.03", "0.1", "0.3"])
    b.xaxis.set_minor_locator(NullLocator())
    b.set_yticks([10, 100, 1000, lag_limit], ["10", "100", "1,000", "3,122.5"])
    b.yaxis.set_minor_locator(NullLocator())
    finite = [row for row in rows if np.isfinite(row["tau_BB_peak_mean"])]
    censored = [row for row in rows if not np.isfinite(row["tau_BB_peak_mean"])]
    b.errorbar([r["k_fb"] for r in finite], [r["tau_BB_peak_mean"] for r in finite],
               yerr=[r["tau_BB_peak_sem"] for r in finite], color=ORANGE, marker="s",
               markersize=5.3, linewidth=1.7, capsize=3, elinewidth=1,
               markeredgecolor="white", markeredgewidth=0.6, label="B density")
    b.errorbar([r["k_fb"] for r in rows], [r["mark_memory_time_mean"] for r in rows],
               yerr=[r["mark_memory_time_sem"] for r in rows], color=INK, marker="D",
               linestyle=(0, (3, 2)), markersize=4.5, linewidth=1.5, capsize=3,
               markeredgecolor="white", markeredgewidth=0.6, label="Site marks")
    b.axhline(lag_limit, color=MUTED, linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
    for row in censored:
        x = row["k_fb"]
        b.plot(x, lag_limit, marker="_", color=ORANGE, markersize=9, markeredgewidth=1.4)
        b.annotate("", xy=(x, 5800), xytext=(x, lag_limit),
                   arrowprops={"arrowstyle": "-|>", "lw": 1.4, "color": ORANGE})
    b.set_xlabel(r"Feedback $k_{\rm fb}$ ($\tau^{-1}$)", labelpad=8)
    b.set_ylabel(r"$1/e$ time ($\tau$)", labelpad=8)
    b.legend(loc="upper left", frameon=False, fontsize=9.5, handlelength=2.1,
             borderaxespad=0.3, labelspacing=0.25)

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
        "per condition. A: each heatmap cell pairs epsilon_BB=1 in its upper-left "
        "and epsilon_BB=1.5 in its lower-right triangle. B: density and site-mark "
        "1/e times at k_off=0.01 and epsilon_BB=1.5; upward arrows are lower "
        "bounds at the 3122.5-tau evaluated lag limit, not finite relaxation estimates. "
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
