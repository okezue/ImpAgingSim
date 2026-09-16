#!/usr/bin/env python3
"""Render the exact nonbonded energy currently implemented in melt/integrator.py.

The tail is zero below the WCA minimum: it has no interior constant extension.
Consequently the reported energy has a jump at r_m.  This figure deliberately
retains that jump rather than drawing a conventional continuous LJ well.

Run from any directory: python scripts/render_potential_figure.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"

INK = "#243345"
MUTED = "#657385"
TEAL = "#007F86"
ORANGE = "#D97943"
PURPLE = "#7556A5"
RULE = "#D9E0E5"

R_MIN = 2.0 ** (1.0 / 6.0)
R_CUT = 2.5


def core(r):
    """Common WCA branch in units sigma=epsilon_core=1."""
    r = np.asarray(r)
    return 4.0 * (r ** -12 - r ** -6) + 1.0


def tail(r, epsilon):
    """Cutoff-shifted attractive branch; call only for r_m <= r <= r_c."""
    r = np.asarray(r)
    shift = 4.0 * epsilon * (R_CUT ** -12 - R_CUT ** -6)
    return 4.0 * epsilon * (r ** -12 - r ** -6) - shift


def bead(ax, xy, label, color):
    ax.add_patch(Circle(xy, 0.125, facecolor=color, edgecolor="white", lw=1.5))
    ax.text(*xy, label, color="white", size=13, weight="bold", ha="center", va="center")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": MUTED,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.linewidth": 0.9,
            "svg.fonttype": "path",
            "svg.hashsalt": "ImpAgingSim-potential-v1",
            "mathtext.fontset": "dejavusans",
            "savefig.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(12, 5), facecolor="white")
    fig.text(0.045, 0.91, "A", size=16, weight="bold")
    fig.text(0.073, 0.91, "Nonbonded pairs", size=15)
    fig.text(0.365, 0.91, "B", size=16, weight="bold")
    fig.text(0.394, 0.91, "Pair energy", size=15)

    pairs = fig.add_axes((0.045, 0.18, 0.265, 0.65))
    pairs.set_aspect("equal")
    pairs.set_xlim(0, 1.05)
    pairs.set_ylim(0, 1.30)
    pairs.axis("off")
    for y, (label1, color1), (label2, color2), strength in (
        (1.08, ("A", TEAL), ("A", TEAL), "1.0"),
        (0.65, ("B", ORANGE), ("B", ORANGE), "1.0"),
        (0.22, ("A", TEAL), ("B", ORANGE), "0.1"),
    ):
        bead(pairs, (0.17, y), label1, color1)
        bead(pairs, (0.50, y), label2, color2)
        pairs.text(0.77, y, rf"$\epsilon={strength}$", va="center", size=13,
                   color=PURPLE if strength == "0.1" else INK)
    fig.text(0.074, 0.12, "Common bead diameter", size=11, color=MUTED)

    ax = fig.add_axes((0.445, 0.19, 0.515, 0.64))
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_color(RULE)
    ax.spines["bottom"].set_color(RULE)
    ax.set_xlim(0.91, 2.72)
    ax.set_ylim(-1.16, 2.10)
    ax.set_xlabel(r"Separation $r/\sigma$", labelpad=11, size=13)
    ax.set_ylabel(r"$U(r)/\epsilon_{\rm core}$", labelpad=11, size=13)
    ax.set_xticks([1.0, 1.5, 2.0, 2.5], ["1.0", "1.5", "2.0", "2.5"])
    ax.set_yticks([-1, 0, 1, 2])
    ax.tick_params(length=3.0, width=0.8, labelsize=11)
    ax.axhline(0, color=RULE, lw=1.0, zorder=0)
    ax.axvline(R_MIN, color=MUTED, lw=0.9, ls=(0, (2.5, 4)), alpha=0.7, zorder=0)
    ax.axvline(R_CUT, color=RULE, lw=0.9, ls=(0, (2.5, 4)), zorder=0)

    r_core = np.linspace(0.91, R_MIN, 350)
    r_tail = np.linspace(R_MIN, R_CUT, 600)
    ax.plot(r_core, core(r_core), color=INK, lw=2.7, solid_capstyle="round", zorder=3)
    ax.plot(r_tail, tail(r_tail, 1.0), color=INK, lw=2.7, solid_capstyle="round", zorder=3)
    ax.plot(r_tail, tail(r_tail, 0.1), color=PURPLE, lw=2.7, solid_capstyle="round", zorder=4)
    ax.plot([R_CUT, 2.72], [0, 0], color=MUTED, lw=1.8, zorder=2)

    # Limiting endpoints make the discontinuity explicit, with no connecting line.
    ax.plot(R_MIN, 0.0, "o", ms=7.0, mfc="white", mec=INK, mew=1.7, zorder=6)
    ax.plot(R_MIN, tail(R_MIN, 1.0), "o", ms=6.5, mfc=INK, mec="white", mew=0.6, zorder=6)
    ax.plot(R_MIN, tail(R_MIN, 0.1), "o", ms=6.5, mfc=PURPLE, mec="white", mew=0.6, zorder=7)

    ax.annotate("Shared repulsive core", xy=(0.990, core(0.990)), xytext=(1.31, 1.59),
                ha="left", va="center", size=12, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8,
                                connectionstyle="angle,angleA=180,angleB=75,rad=0"))
    ax.text(1.58, 0.68, "Pair-specific attraction", size=12, color=MUTED)
    ax.annotate(r"$\mathrm{AB}$", xy=(1.71, tail(1.71, 0.1)), xytext=(1.91, 0.26),
                color=PURPLE, size=13, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=PURPLE, lw=0.85))
    ax.annotate(r"$\mathrm{AA}=\mathrm{BB}$", xy=(1.50, tail(1.50, 1.0)), xytext=(1.73, -0.65),
                color=INK, size=13, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.85))
    ax.text(R_MIN + 0.035, 1.97, r"$r_m=2^{1/6}\sigma$", size=11, color=MUTED, va="top")
    ax.text(R_CUT - 0.035, 1.97, r"$r_c=2.5\sigma$", size=11, color=MUTED,
            ha="right", va="top")
    ax.text(1.13, -1.075, "Tail onset", size=10, color=MUTED, ha="left", va="center")

    # Fixed canvas preserves a consistent manuscript-style aspect ratio in the README.
    description = (
        "Exact implemented split nonbonded energy. Universal WCA core below "
        "2^(1/6) sigma; pair-specific LJ attraction above that radius, shifted "
        "to zero at 2.5 sigma. The step-gated tail produces an energy jump at "
        "its onset. Neighboring bonded beads are excluded."
    )
    fig.savefig(OUT / "interaction-potential.svg", metadata={"Date": None, "Title": "Split nonbonded interactions", "Description": description})
    fig.savefig(OUT / "interaction-potential.png", dpi=180, metadata={"Description": description})
    plt.close(fig)
    print(OUT / "interaction-potential.svg")
    print(OUT / "interaction-potential.png")


if __name__ == "__main__":
    main()
