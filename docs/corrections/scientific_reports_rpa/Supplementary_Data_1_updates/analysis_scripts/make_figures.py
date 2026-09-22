#!/usr/bin/env python3
"""Generate the manuscript-facing main and supplementary static figures.

Run ``rebuild_data.py`` first.  Every plotted value is then read from a CSV in
``source_data``; the plotting layer performs no hidden scientific calculation.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path
from uuid import uuid4

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import linregress


HERE = Path(__file__).resolve().parent
SD = HERE / "source_data"
PDF = HERE / "figures" / "pdf"
PNG = HERE / "figures" / "png"
SUBMISSION_PDF = HERE / "figures" / "submission" / "pdf"
SUBMISSION_PNG = HERE / "figures" / "submission" / "png"
PDF.mkdir(parents=True, exist_ok=True)
PNG.mkdir(parents=True, exist_ok=True)
SUBMISSION_PDF.mkdir(parents=True, exist_ok=True)
SUBMISSION_PNG.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
PURPLE = "#7A5195"
SKY = "#56B4E9"
MAGENTA = "#CC79A7"
GRAY = "#5C6770"
LIGHTGRAY = "#C7CDD1"
BLACK = "#1B1F23"

KCOL = {0.0: BLUE, 0.5: ORANGE, 0.6: ORANGE, 1.0: PURPLE}
COND_COL = {
    "baseline": BLUE,
    "high-density": GREEN,
    "soft tail/T": PURPLE,
    "short-chain": RED,
}
SEQ_COL = {
    "correlated": GREEN,
    "random": GRAY,
    "block": RED,
    "alternating": SKY,
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.4,
        "axes.labelsize": 8.6,
        "axes.titlesize": 9.0,
        "axes.titlepad": 8.0,
        "xtick.labelsize": 7.6,
        "ytick.labelsize": 7.6,
        "legend.fontsize": 7.4,
        "axes.linewidth": 1.0,
        "axes.facecolor": "white",
        "axes.edgecolor": BLACK,
        "grid.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "lines.linewidth": 1.45,
        "lines.markeredgewidth": 1.0,
        "patch.linewidth": 1.0,
        "hatch.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "savefig.transparent": False,
        "figure.facecolor": "white",
        "figure.edgecolor": "white",
    }
)


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(SD / name)


def panel(ax: plt.Axes, label: str) -> None:
    label = label.lower()
    ax.annotate(
        label,
        xy=(0, 1),
        xycoords="axes fraction",
        xytext=(-10, 20),
        textcoords="offset points",
        ha="right",
        va="bottom",
        fontsize=11.5,
        fontweight="bold",
        clip_on=False,
    )


def polish(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis="y", color="#D7DCE0", lw=1.0, alpha=0.65, zorder=0)
    ax.set_axisbelow(True)


def raw_points_at(
    ax: plt.Axes,
    x: float,
    values,
    *,
    color: str,
    jitter: float = 0.04,
    log_x: bool = False,
    alpha: float = 0.32,
    zorder: float = 3.2,
) -> None:
    """Overlay a restrained deterministic seed cloud without adding thin edges."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return
    offsets = np.zeros(1) if len(values) == 1 else np.linspace(-0.5, 0.5, len(values))
    if log_x:
        x_values = float(x) * np.exp(offsets * jitter)
    else:
        x_values = float(x) + offsets * jitter
    ax.scatter(
        x_values,
        values,
        s=11,
        facecolor=color,
        edgecolors="none",
        linewidths=0,
        alpha=alpha,
        zorder=zorder,
    )


def overlay_seed_points(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    *,
    color: str,
    jitter: float = 0.04,
    log_x: bool = False,
) -> None:
    """Overlay one point per independent seed at every plotted condition."""
    clean = data[[x_column, y_column] + (["seed"] if "seed" in data else [])].dropna()
    for x_value, group in clean.groupby(x_column, sort=True):
        if "seed" in group:
            group = group.sort_values("seed")
        raw_points_at(
            ax,
            float(x_value),
            group[y_column].to_numpy(float),
            color=color,
            jitter=jitter,
            log_x=log_x,
        )


def publication_colorbar(fig: plt.Figure, mappable, **kwargs):
    """Keep colour scales vector and enforce the journal's one-point rule."""
    cb = fig.colorbar(mappable, **kwargs)
    if cb.solids is not None:
        cb.solids.set_rasterized(False)
    cb.outline.set_linewidth(1.0)
    cb.ax.tick_params(which="both", width=1.0)
    for spine in cb.ax.spines.values():
        spine.set_linewidth(1.0)
    return cb


def legend_above(
    ax: plt.Axes,
    *,
    ncol: int = 1,
    anchor: tuple[float, float] = (0.5, 1.015),
    **kwargs,
):
    """Place a legend outside the plotting field so it cannot hide data."""
    return ax.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=anchor,
        borderaxespad=0,
        ncol=ncol,
        columnspacing=1.0,
        handletextpad=0.45,
        **kwargs,
    )


def finish(fig: plt.Figure, name: str) -> None:
    pdf_path = PDF / f"{name}.pdf"
    png_path = PNG / f"{name}.png"
    write_id = uuid4().hex
    pdf_tmp = PDF / f".{name}.{write_id}.pdf.tmp"
    png_tmp = PNG / f".{name}.{write_id}.png.tmp"
    try:
        fig.patch.set_facecolor("white")
        fig.patch.set_alpha(1.0)
        for ax in fig.axes:
            ax.set_facecolor("white")
        fig.savefig(
            pdf_tmp,
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.04,
            facecolor="white",
            edgecolor="white",
            transparent=False,
        )
        fig.savefig(
            png_tmp,
            format="png",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.04,
            facecolor="white",
            edgecolor="white",
            transparent=False,
        )
        pdf_tmp.replace(pdf_path)
        png_tmp.replace(png_path)
    finally:
        pdf_tmp.unlink(missing_ok=True)
        png_tmp.unlink(missing_ok=True)
        plt.close(fig)


def write_submission_aliases() -> None:
    """Copy the generated figures to manuscript-facing stable filenames."""
    aliases = {
        "Fig1_sequence_construction": "Fig1_sequence_statistics",
        "Fig2_kappa_controls": "Fig2_condition_controls",
        "Fig3_pi_kappa_response": "Fig3_pi_kappa_response",
        "Fig4_sequence_response_theory": "Fig4_finite_predictor_rpa",
        "Fig5_cross_attraction": "Fig5_cross_attraction",
        "FigS1_sequence_validation": "FigS2_sequence_validation",
        "FigS2_form_factor_limits": "FigS7_rpa_form_factors",
        "FigS3_peak_binning": "Fig6_peak_binning_sensitivity",
        "FigS4_cic_sensitivity": "FigS1_representative_spectra",
        "FigS5_kappa_controls": "FigS6_condition_diagnostics",
        "FigS6_composition_confounding": "FigS4_composition_observable_limitation",
        "FigS7_interaction_statistics": "FigS5_interaction_controls",
    }
    for output_dir, source_dir, suffix in (
        (SUBMISSION_PDF, PDF, ".pdf"),
        (SUBMISSION_PNG, PNG, ".png"),
    ):
        for stale in output_dir.glob(f"*{suffix}"):
            stale.unlink()
        for destination, source in aliases.items():
            destination_path = output_dir / f"{destination}{suffix}"
            temporary_path = output_dir / f".{destination}.{uuid4().hex}{suffix}.tmp"
            try:
                shutil.copy2(source_dir / f"{source}{suffix}", temporary_path)
                temporary_path.replace(destination_path)
            finally:
                temporary_path.unlink(missing_ok=True)


def sync_journal_destinations() -> None:
    """Synchronize stable publication figures without touching Supplementary Fig. S8."""
    package_root = HERE.parents[2]
    journal = package_root / "journal_upload"
    main_names = [
        "Fig1_sequence_construction",
        "Fig2_kappa_controls",
        "Fig3_pi_kappa_response",
        "Fig4_sequence_response_theory",
        "Fig5_cross_attraction",
    ]
    supplementary_names = [
        "FigS1_sequence_validation",
        "FigS2_form_factor_limits",
        "FigS3_peak_binning",
        "FigS4_cic_sensitivity",
        "FigS5_kappa_controls",
        "FigS6_composition_confounding",
        "FigS7_interaction_statistics",
    ]
    destinations = [
        (SUBMISSION_PDF, journal / "figures" / "main_pdf", main_names, ".pdf"),
        (SUBMISSION_PNG, journal / "figures" / "main_png", main_names, ".png"),
        (SUBMISSION_PDF, journal / "manuscript" / "figures", main_names, ".pdf"),
        (SUBMISSION_PDF, journal / "manuscript_source" / "figures", main_names, ".pdf"),
        (
            SUBMISSION_PDF,
            journal / "supplementary" / "figures",
            supplementary_names,
            ".pdf",
        ),
    ]
    for source_dir, destination_dir, names, suffix in destinations:
        destination_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            source = source_dir / f"{name}{suffix}"
            destination = destination_dir / source.name
            temporary = destination_dir / f".{name}.{uuid4().hex}{suffix}.tmp"
            try:
                shutil.copy2(source, temporary)
                temporary.replace(destination)
            finally:
                temporary.unlink(missing_ok=True)


def clean_figure_outputs() -> None:
    """Remove only figures generated by this script before redrawing."""
    for output_dir, suffix in (
        (PDF, ".pdf"),
        (PNG, ".png"),
        (SUBMISSION_PDF, ".pdf"),
        (SUBMISSION_PNG, ".png"),
    ):
        for path in output_dir.glob(f"*{suffix}"):
            path.unlink()
        for path in output_dir.glob(".*.tmp"):
            path.unlink()


def matrix(df: pd.DataFrame, value: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    piv = df.pivot(index="pi", columns="kappa", values=value).sort_index().sort_index(axis=1)
    return piv.columns.to_numpy(float), piv.index.to_numpy(float), piv.to_numpy(float)


def centers_to_edges(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    mid = 0.5 * (v[:-1] + v[1:])
    first = v[0] - (mid[0] - v[0])
    last = v[-1] + (v[-1] - mid[-1])
    return np.r_[first, mid, last]


def categorical_y_heatmap(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    **kwargs,
):
    """Draw sampled pi rows without implying equal numeric spacing in pi."""
    x_edges = centers_to_edges(x)
    y_edges = np.arange(len(y) + 1, dtype=float) - 0.5
    image = ax.pcolormesh(x_edges, y_edges, z, shading="flat", **kwargs)
    ax.set_yticks(np.arange(len(y)), [f"{value:g}" for value in y])
    ax.set_ylim(-0.5, len(y) - 0.5)
    return image


def fig1_sequence_statistics() -> None:
    strips = read("Fig1_sequence_strips.csv")
    cov = read("Fig1_sequence_covariance.csv")
    bound = read("Fig1_chain_boundary_check.csv")
    form = read("Fig1_finite_chain_form_factor_q0.csv")

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.25), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    cov_plot = (
        cov.groupby(["kappa", "lag"], as_index=False)
        .agg(
            covariance_mean=("covariance_seed_mean", "mean"),
            covariance_sem=("covariance_seed_mean", "sem"),
            covariance_theory=("covariance_theory", "first"),
            n_seeds=("seed", "size"),
        )
    )
    bound_plot = (
        bound.groupby(["kappa", "comparison"], as_index=False)
        .agg(
            covariance_mean=("covariance_seed_mean", "mean"),
            covariance_sem=("covariance_seed_mean", "sem"),
            theory=("theory", "first"),
            n_seeds=("seed", "size"),
        )
    )

    # A: exact seed-1 reconstruction of independent static-scan chain draws.
    row_labels = []
    image_rows = []
    for kappa in (0.0, 0.6, 1.0):
        for chain in range(4):
            sub = strips[(strips.kappa == kappa) & (strips.chain_index == chain)].sort_values(
                "bead_index"
            )
            image_rows.append(sub.type_numeric.to_numpy())
            row_labels.append(rf"$\kappa={kappa:g}$, c{chain + 1}")
    cmap = ListedColormap([BLUE, RED])
    image_rows = np.asarray(image_rows)
    ax_a.pcolormesh(
        np.arange(image_rows.shape[1] + 1) + 0.5,
        np.arange(image_rows.shape[0] + 1) - 0.5,
        image_rows,
        shading="flat",
        cmap=cmap,
        vmin=0,
        vmax=1,
        rasterized=False,
    )
    for y in (3.5, 7.5):
        ax_a.axhline(y, color="white", lw=2.5)
    ax_a.set_yticks(np.arange(len(row_labels)), row_labels)
    ax_a.set_xlim(0.5, 40.5)
    ax_a.set_xticks([1, 10, 20, 30, 40], [1, 10, 20, 30, 40])
    ax_a.set_xlabel("monomer index")
    legend_above(
        ax_a,
        handles=[Patch(color=RED, label="A"), Patch(color=BLUE, label="B")],
        ncol=2,
    )
    panel(ax_a, "A")

    # B: within-chain covariance, actual chains versus exact two-stage theory.
    for kappa in (0.0, 0.6, 1.0):
        sub = cov_plot[cov_plot.kappa == kappa].sort_values("lag")
        raw = cov[cov.kappa == kappa].sort_values(["lag", "seed"])
        color = KCOL[kappa]
        overlay_seed_points(ax_b, raw, "lag", "covariance_seed_mean", color=color, jitter=0.18)
        ax_b.errorbar(
            sub.lag,
            sub.covariance_mean,
            yerr=sub.covariance_sem,
            fmt="o",
            ms=3.0,
            color=color,
            capsize=1.6,
            label=rf"$\kappa={kappa:g}$",
            zorder=4,
        )
        ax_b.plot(sub.lag, sub.covariance_theory, color=color, lw=1.35)
    ax_b.axhline(0, color=BLACK, lw=1.0)
    ax_b.set_xlabel(r"within-chain lag $\ell$")
    ax_b.set_ylabel(r"$\langle\delta\psi_i\delta\psi_{i+\ell}\rangle$")
    legend_above(ax_b, ncol=3)
    polish(ax_b)
    panel(ax_b, "B")

    # C: within-chain versus chain-boundary covariance.
    x = np.arange(3)
    width = 0.34
    for j, comparison in enumerate(("within-chain lag 1", "across-chain boundary")):
        sub = bound_plot[bound_plot.comparison == comparison].sort_values("kappa")
        ax_c.bar(
            x + (j - 0.5) * width,
            sub.covariance_mean,
            width,
            yerr=sub.covariance_sem,
            color=GREEN if j == 0 else LIGHTGRAY,
            edgecolor=BLACK,
            linewidth=1.0,
            capsize=2.2,
            label=comparison,
            zorder=3,
        )
        ax_c.scatter(
            x + (j - 0.5) * width,
            sub.theory,
            marker="_",
            s=80,
            color=BLACK,
            zorder=4,
        )
        raw = bound[bound.comparison == comparison]
        for i, kappa in enumerate((0.0, 0.6, 1.0)):
            values = raw[np.isclose(raw.kappa, kappa)].sort_values("seed").covariance_seed_mean
            raw_points_at(
                ax_c,
                x[i] + (j - 0.5) * width,
                values,
                color=BLACK,
                jitter=0.08,
                alpha=0.38,
                zorder=4.2,
            )
    ax_c.axhline(0, color=BLACK, lw=1.0)
    ax_c.set_xticks(x, [r"$0$", r"$0.6$", r"$1$"])
    ax_c.set_xlabel(r"$\kappa$")
    ax_c.set_ylabel("sequence covariance")
    legend_above(ax_c, ncol=2)
    polish(ax_c)
    panel(ax_c, "C")

    # D: independent amplitude and persistence control.
    for pi, color in zip((0.70, 0.85, 0.95, 0.99), (SKY, GREEN, ORANGE, PURPLE)):
        sub = form[np.isclose(form.pi, pi)]
        ax_d.plot(sub.kappa, sub.P_N_q0_normalized, color=color, label=rf"$\pi={pi:g}$")
    ax_d.set_yscale("log")
    ax_d.set_xlabel(r"backbone-expression probability $\kappa$")
    ax_d.set_ylabel(r"finite-chain $P_N(0)$")
    legend_above(ax_d, ncol=2)
    polish(ax_d)
    panel(ax_d, "D")
    finish(fig, "Fig1_sequence_statistics")


def fig2_condition_controls() -> None:
    df = read("FigS6_condition_control_summary.csv").copy()
    runs = read("FigS6_condition_control_runs.csv").copy()
    df["condition"] = df.condition.replace(
        {"dense": "high-density", "soft": "soft tail/T", "short chain": "short-chain"}
    )
    runs["condition"] = runs.condition.replace(
        {"dense": "high-density", "soft": "soft tail/T", "short chain": "short-chain"}
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.05), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    for condition in ("baseline", "high-density", "soft tail/T", "short-chain"):
        sub = df[df.condition == condition].sort_values("kappa")
        raw = runs[runs.condition == condition].sort_values(["kappa", "seed"])
        color = COND_COL[condition]
        overlay_seed_points(ax_a, raw, "kappa", "C_archive_N", color=color, jitter=0.018)
        ax_a.errorbar(
            sub.kappa,
            sub.C_archive_N_mean,
            yerr=sub.C_archive_N_sem,
            marker="o",
            ms=3.6,
            capsize=2,
            color=color,
            label=condition,
            zorder=4,
        )
        c0 = float(sub.loc[sub.kappa.idxmin(), "C_archive_N_mean"])
        raw = raw.assign(C_relative=raw.C_archive_N / c0)
        overlay_seed_points(ax_b, raw, "kappa", "C_relative", color=color, jitter=0.018)
        ax_b.errorbar(
            sub.kappa,
            sub.C_archive_N_mean / c0,
            yerr=sub.C_archive_N_sem / c0,
            marker="o",
            ms=3.6,
            capsize=2,
            color=color,
            zorder=4,
        )
        ax_c.plot(sub.kappa, sub.k_star_mode, marker="o", ms=3.6, color=color)
        overlay_seed_points(ax_d, raw, "kappa", "Rg", color=color, jitter=0.018)
        ax_d.errorbar(
            sub.kappa,
            sub.Rg_mean,
            yerr=sub.Rg_sem,
            marker="o",
            ms=3.6,
            capsize=2,
            color=color,
            zorder=4,
        )
    ax_a.set_yscale("log")
    ax_a.set_ylabel(r"spectral peak amplitude $\widehat C_N$")
    ax_b.set_ylabel(r"$\widehat C_N(\kappa)/\widehat C_N(0)$")
    ax_c.set_ylabel(r"modal $k^*\sigma$")
    ax_d.set_ylabel(r"mean chain $R_g/\sigma$")
    for ax in axes.ravel():
        ax.set_xlabel(r"$\kappa$")
        polish(ax)
    legend_above(ax_a, ncol=2)
    for ax, lab in zip(axes.ravel(), "ABCD"):
        panel(ax, lab)
    finish(fig, "Fig2_condition_controls")


def fig3_pi_kappa_response() -> None:
    df = read("Fig2_pi_kappa_summary.csv")
    runs = read("Fig2_pi_kappa_runs.csv")
    x, y, zc = matrix(df, "C_archive_N_mean")
    _, _, zk = matrix(df, "k_star_mode")
    _, _, zf = matrix(df, "first_mixed_bin_fraction")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.25), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    im = categorical_y_heatmap(
        ax_a,
        x,
        y,
        zc,
        cmap="viridis",
        norm=LogNorm(zc.min(), zc.max()),
    )
    cb = publication_colorbar(fig, im, ax=ax_a, pad=0.02)
    cb.set_label(r"mean $\widehat C_N$")
    ax_a.set_xlabel(r"$\kappa$")
    ax_a.set_ylabel(r"sampled $\pi$")
    panel(ax_a, "A")

    levels = sorted(np.unique(zk))
    bounds = np.r_[np.array(levels) - 0.08, levels[-1] + 0.08]
    cmap_k = ListedColormap([SKY, GREEN, ORANGE, PURPLE])
    im2 = categorical_y_heatmap(
        ax_b,
        x,
        y,
        zk,
        cmap=cmap_k,
        norm=BoundaryNorm(bounds, cmap_k.N),
    )
    cb2 = publication_colorbar(fig, im2, ax=ax_b, pad=0.02, ticks=levels)
    cb2.set_label(r"modal $k^*\sigma$")
    ax_b.set_xlabel(r"$\kappa$")
    ax_b.set_ylabel(r"sampled $\pi$")
    panel(ax_b, "B")

    for pi, color in ((0.50, BLUE), (0.85, ORANGE), (0.99, PURPLE)):
        sub = df[np.isclose(df.pi, pi)].sort_values("kappa")
        raw = runs[np.isclose(runs.pi, pi)].sort_values(["kappa", "seed"])
        overlay_seed_points(ax_c, raw, "kappa", "C_archive_N", color=color, jitter=0.018)
        ax_c.errorbar(
            sub.kappa,
            sub.C_archive_N_mean,
            yerr=sub.C_archive_N_sem,
            marker="o",
            ms=3.8,
            capsize=2,
            color=color,
            label=rf"$\pi={pi:g}$",
            zorder=4,
        )
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"$\kappa$")
    ax_c.set_ylabel(r"$\widehat C_N$")
    legend_above(ax_c, ncol=3)
    polish(ax_c)
    panel(ax_c, "C")

    im4 = categorical_y_heatmap(ax_d, x, y, zf, cmap="magma", vmin=0, vmax=1)
    cb4 = publication_colorbar(fig, im4, ax=ax_d, pad=0.02)
    cb4.set_label(r"fraction at $k^*\sigma=0.4284$")
    ax_d.set_xlabel(r"$\kappa$")
    ax_d.set_ylabel(r"sampled $\pi$")
    panel(ax_d, "D")
    finish(fig, "Fig3_pi_kappa_response")


def fig4_finite_predictor_rpa() -> None:
    pts = read("Fig6_predictor_points.csv")
    reg = read("Fig6_predictor_regression_audit.csv")
    form = read("Fig6_finite_chain_form_factor.csv")
    rpa = read("Fig6_rpa_q0_diagnostic.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.1), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    primary = pts[pts.primary_regression_included.astype(bool)]
    sc = ax_a.scatter(
        primary.inverse_P_N_kstar,
        primary.inverse_C_archive_N,
        c=primary.pi,
        cmap="plasma",
        s=25,
        edgecolor="white",
        linewidth=1.0,
        zorder=3,
    )
    lr = linregress(primary.inverse_P_N_kstar, primary.inverse_C_archive_N)
    xx = np.linspace(primary.inverse_P_N_kstar.min(), primary.inverse_P_N_kstar.max(), 100)
    ax_a.plot(xx, lr.intercept + lr.slope * xx, color=BLACK, lw=1.2)
    publication_colorbar(fig, sc, ax=ax_a, pad=0.02, label=r"$\pi$")
    ax_a.set_xlabel(r"$1/P_N(k^*)$")
    ax_a.set_ylabel(r"$1/\langle\widehat C_N\rangle$")
    polish(ax_a)
    panel(ax_a, "A")

    arch = reg[reg.dataset.str.startswith("available")]
    xpos = np.arange(2)
    vals = [
        float(arch[arch.predictor == "finite q=k_star"].r_squared.iloc[0]),
        float(arch[arch.predictor == "q=0"].r_squared.iloc[0]),
    ]
    bars = ax_b.bar(xpos, vals, color=[GREEN, LIGHTGRAY], edgecolor=BLACK, linewidth=1.0, zorder=3)
    ax_b.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)
    ax_b.set_xticks(xpos, [r"finite $q=k^*$", r"$q=0$"])
    ax_b.set_ylim(0, 1.05)
    ax_b.set_ylabel(r"$R^2$")
    polish(ax_b)
    panel(ax_b, "B")

    for pi, color in ((0.50, BLUE), (0.85, ORANGE), (0.99, PURPLE)):
        sub = form[(np.isclose(form.pi, pi)) & (np.isclose(form.kappa, 1.0))]
        ax_c.plot(sub.q_sigma, sub.P_N_q_normalized, color=color, label=rf"$\pi={pi:g}$")
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"$q\sigma$")
    ax_c.set_ylabel(r"$P_N(q)$")
    legend_above(ax_c, ncol=3)
    polish(ax_c)
    panel(ax_c, "C")

    ideal = float(rpa[rpa.quantity.str.startswith("ideal")].value.iloc[0])
    interaction = -float(rpa[rpa.quantity.str.startswith("interaction")].value.iloc[0])
    total = float(rpa[rpa.quantity.isin([
        "bare composition inverse stiffness", "bare scalar RPA denominator"
    ])].value.iloc[0])
    rbars = ax_d.bar(
        np.arange(3),
        [ideal, interaction, total],
        color=[BLUE, RED, PURPLE],
        edgecolor=BLACK,
        linewidth=1.0,
        zorder=3,
    )
    ax_d.bar_label(rbars, fmt="%.2f", padding=3, fontsize=8)
    ax_d.axhline(0, color=BLACK, lw=1.0)
    ax_d.set_xticks(np.arange(3), ["ideal", r"$-\chi_{\rm bare}/2$", "bare total"])
    ax_d.set_ylim(-5.55, 1.35)
    ax_d.set_ylabel(r"bare composition stiffness ($q\to0$)")
    polish(ax_d)
    panel(ax_d, "D")
    finish(fig, "Fig4_finite_predictor_rpa")


def fig5_cross_attraction() -> None:
    runs = read("Fig5_interaction_runs.csv")
    summ = read("Fig5_interaction_summary.csv")
    tests = read("Fig5_interaction_tests.csv")
    corr = summ[summ.sequence == "correlated"].sort_values("eps_AB")
    cruns = runs[runs.sequence == "correlated"].copy()
    permutation = tests[
        (tests.analysis == "available processed archived estimator")
        & tests.test.str.startswith("1,000,000", na=False)
    ].iloc[0]
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.1), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    overlay_seed_points(
        ax_a,
        cruns,
        "eps_AB",
        "C_archive_N",
        color=GREEN,
        jitter=0.08,
        log_x=True,
    )
    ax_a.errorbar(
        corr.eps_AB,
        corr.C_archive_N_mean,
        yerr=corr.C_archive_N_sem,
        marker="o",
        ms=3.8,
        capsize=2,
        color=GREEN,
        zorder=4,
    )
    ax_a.set_xscale("log")
    ax_a.set_xlabel(r"cross attraction $\epsilon_{AB}/\epsilon_{AA}$")
    ax_a.set_ylabel(r"$\widehat C_N$")
    polish(ax_a)
    panel(ax_a, "A")

    low = cruns[cruns.eps_AB <= 0.2]
    for seed, sub in low.groupby("seed"):
        sub = sub.sort_values("eps_AB")
        ax_b.plot(sub.eps_AB, sub.C_archive_N, color=LIGHTGRAY, lw=1.0, alpha=0.9)
    low_s = corr[corr.eps_AB <= 0.2]
    ax_b.errorbar(
        low_s.eps_AB,
        low_s.C_archive_N_mean,
        yerr=low_s.C_archive_N_sem,
        marker="o",
        ms=3.8,
        capsize=2,
        color=GREEN,
        zorder=4,
    )
    ax_b.set_xscale("log")
    ax_b.set_xlabel(r"$\epsilon_{AB}/\epsilon_{AA}$")
    ax_b.set_ylabel(r"$\widehat C_N$")
    polish(ax_b)
    panel(ax_b, "B")

    paired = cruns[cruns.eps_AB.isin([0.1, 0.8])].pivot(index="seed", columns="eps_AB", values="C_archive_N")
    for seed, row in paired.iterrows():
        ax_c.plot([0, 1], [row[0.1], row[0.8]], marker="o", ms=3.5, color=GRAY, alpha=0.75)
    endpoint = corr[corr.eps_AB.isin([0.1, 0.8])].sort_values("eps_AB")
    ax_c.errorbar(
        [0, 1],
        endpoint.C_archive_N_mean,
        yerr=endpoint.C_archive_N_sem,
        color=GREEN,
        marker="D",
        ms=5,
        lw=2.0,
        capsize=3,
        zorder=5,
    )
    ax_c.set_xticks([0, 1], [r"$\epsilon_{AB}=0.1$", r"$\epsilon_{AB}=0.8$"])
    ax_c.set_ylabel(r"$\widehat C_N$")
    polish(ax_c)
    panel(ax_c, "C")

    use = summ[summ.eps_AB.isin([0.1, 0.8])].copy()
    seq_order = ["alternating", "random", "block", "correlated"]
    x = np.arange(len(seq_order))
    width = 0.34
    for j, eps in enumerate((0.1, 0.8)):
        sub = use[np.isclose(use.eps_AB, eps)].set_index("sequence").loc[seq_order]
        ax_d.bar(
            x + (j - 0.5) * width,
            sub.C_archive_N_mean,
            width,
            yerr=sub.C_archive_N_sem,
            capsize=2,
            color=SKY if eps == 0.1 else ORANGE,
            edgecolor=BLACK,
            linewidth=1.0,
            label=rf"$\epsilon_{{AB}}={eps:g}$",
            zorder=3,
        )
        for i, sequence in enumerate(seq_order):
            values = runs[
                (runs.sequence == sequence) & np.isclose(runs.eps_AB, eps)
            ].sort_values("seed").C_archive_N
            raw_points_at(
                ax_d,
                x[i] + (j - 0.5) * width,
                values,
                color=BLACK,
                jitter=0.08,
                alpha=0.38,
                zorder=4.2,
            )
    ax_d.set_yscale("log")
    ax_d.set_xticks(x, seq_order, rotation=18, ha="right")
    ax_d.set_ylabel(r"$\widehat C_N$")
    legend_above(ax_d, ncol=2)
    polish(ax_d)
    panel(ax_d, "D")
    finish(fig, "Fig5_cross_attraction")


def fig6_peak_binning_sensitivity() -> None:
    fig, ax = plt.subplots(figsize=(7.0, 3.0), constrained_layout=True)
    items = [
        ("symmetric peak rule", 4, PURPLE),
        ("CIC peak rule", 1, GREEN),
        ("symmetric first-bin calls", 98, ORANGE),
        (r"$C_A$ first-bin calls", 99, BLUE),
    ]
    vals = [x[1] for x in items]
    bars = ax.bar(np.arange(len(items)), vals, color=[x[2] for x in items], edgecolor=BLACK, linewidth=1.0, zorder=3)
    ax.bar_label(bars, labels=[f"{v}/240" for v in vals], padding=3, fontsize=8)
    ax.set_xticks(np.arange(len(items)), [x[0] for x in items], rotation=10, ha="right")
    ax.set_ylabel("affected or selected runs")
    polish(ax)
    finish(fig, "Fig6_peak_binning_sensitivity")


def figs1_representative_spectra() -> None:
    s = read("Fig3_representative_spectra_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.9), constrained_layout=True)
    ax_a, ax_b, ax_c = axes.ravel()

    for kappa in (0.0, 0.5, 1.0):
        sub = s[s.kappa == kappa].sort_values("q_sigma")
        color = KCOL[kappa]
        ax_a.plot(sub.q_sigma, sub.C_N_direct_mean, color=color, marker="o", ms=2.2, label=rf"$\kappa={kappa:g}$")
        lo = np.clip(sub.C_N_direct_mean - sub.C_N_direct_sd_over_frames, 1e-5, None)
        hi = sub.C_N_direct_mean + sub.C_N_direct_sd_over_frames
        ax_a.fill_between(sub.q_sigma, lo, hi, color=color, alpha=0.12, linewidth=0)
        ax_a.plot(sub.q_sigma, sub.C_N_CIC_mean, color=color, ls="--", lw=1.0)
    ax_a.plot([], [], color=GRAY, ls="--", lw=1.0, label="raw CIC")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"$q\sigma$")
    ax_a.set_ylabel(r"$C_N(q)$")
    legend_above(ax_a, ncol=2)
    polish(ax_a)

    sub = s[s.kappa == 1.0].sort_values("q_sigma")
    for col, label, color in (
        ("S_AA_N_direct_mean", r"$S_{AA}^{(N)}$", RED),
        ("S_BB_N_direct_mean", r"$S_{BB}^{(N)}$", BLUE),
        ("S_AB_N_direct_mean", r"$S_{AB}^{(N)}$", GREEN),
    ):
        ax_b.plot(sub.q_sigma, sub[col], marker="o", ms=2.2, color=color, label=label)
    ax_b.set_yscale("symlog", linthresh=0.5)
    ax_b.axhline(0, color=BLACK, lw=1.0)
    ax_b.set_xlabel(r"exact shell $q\sigma$")
    ax_b.set_ylabel("number-normalized partial spectrum")
    legend_above(ax_b, ncol=3)
    polish(ax_b)

    for kappa in (0.0, 0.5, 1.0):
        sub = s[s.kappa == kappa].sort_values("q_sigma")
        ax_c.plot(sub.q_sigma, sub.CIC_raw_bias_percent, color=KCOL[kappa], label=rf"$\kappa={kappa:g}$")
    ax_c.axhline(0, color=BLACK, lw=1.0)
    ax_c.set_xlabel(r"$q\sigma$")
    ax_c.set_ylabel("CIC bias relative to direct (%)")
    legend_above(ax_c, ncol=3)
    polish(ax_c)
    for ax, lab in zip(axes.ravel(), "ABC"):
        panel(ax, lab)
    finish(fig, "FigS1_representative_spectra")


def figs2_sequence_validation() -> None:
    cov = read("Fig1_sequence_covariance.csv")
    strips = read("Fig1_sequence_strips.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9), constrained_layout=True)
    ax_a, ax_b = axes.ravel()
    cov_plot = cov.groupby(["kappa", "lag"], as_index=False).agg(observed=("covariance_seed_mean", "mean"), se=("covariance_seed_mean", "sem"), theory=("covariance_theory", "first"))
    for kappa in (0.0, 0.6, 1.0):
        sub = cov_plot[cov_plot.kappa == kappa]
        raw = cov[cov.kappa == kappa].copy()
        raw["residual"] = raw.covariance_seed_mean - raw.covariance_theory
        overlay_seed_points(ax_a, raw, "lag", "residual", color=KCOL[kappa], jitter=0.18)
        ax_a.errorbar(sub.lag, sub.observed - sub.theory, yerr=sub.se, marker="o", ms=2.5, capsize=1.5, color=KCOL[kappa], label=rf"$\kappa={kappa:g}$", zorder=4)
    ax_a.axhline(0, color=BLACK, lw=1.0)
    ax_a.set_xlabel("lag")
    ax_a.set_ylabel("observed minus theory")
    legend_above(ax_a, ncol=3)
    polish(ax_a)

    realized = strips.groupby("kappa").type_numeric.mean().reset_index()
    ax_b.bar(realized.kappa.astype(str), realized.type_numeric, color=[BLUE, ORANGE, PURPLE], edgecolor=BLACK, linewidth=1.0)
    ax_b.axhline(0.5, color=BLACK, ls=":", lw=1.0)
    ax_b.set_ylim(0, 1.05)
    ax_b.set_xlabel(r"$\kappa$")
    ax_b.set_ylabel("A fraction in 4 displayed chains")
    polish(ax_b)
    for ax, lab in zip(axes.ravel(), "AB"):
        panel(ax, lab)
    finish(fig, "FigS2_sequence_validation")


def figs3_pi_kappa_replicates() -> None:
    runs = read("Fig2_pi_kappa_runs.csv")
    summ = read("Fig2_pi_kappa_summary.csv")
    x, y, zsem = matrix(summ, "C_archive_N_sem")
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.9), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    rng = np.random.default_rng(123)
    for pi, color in ((0.50, BLUE), (0.85, ORANGE), (0.99, PURPLE)):
        sub = runs[np.isclose(runs.pi, pi)]
        ax_a.scatter(sub.kappa + rng.normal(0, 0.008, len(sub)), sub.C_archive_N, s=10, alpha=0.5, color=color, label=rf"$\pi={pi:g}$")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"$\kappa$")
    ax_a.set_ylabel(r"per-run $\widehat C_N$")
    legend_above(ax_a, ncol=3)
    polish(ax_a)

    im = categorical_y_heatmap(ax_b, x, y, zsem, cmap="magma")
    publication_colorbar(fig, im, ax=ax_b, pad=0.02, label="seed SEM")
    ax_b.set_xlabel(r"$\kappa$")
    ax_b.set_ylabel(r"sampled $\pi$")

    qvals = sorted(runs.k_star.unique())
    agg = runs.groupby(["pi", "k_star"]).size().unstack(fill_value=0).reindex(columns=qvals, fill_value=0)
    bottom = np.zeros(len(agg))
    for q, color in zip(qvals, (SKY, GREEN, ORANGE, PURPLE)):
        vals = agg[q].to_numpy()
        ax_c.bar(agg.index.astype(str), vals, bottom=bottom, color=color, label=f"{q:.3f}")
        bottom += vals
    ax_c.set_xlabel(r"$\pi$")
    ax_c.set_ylabel(r"runs across all $\kappa$")
    legend_above(ax_c, ncol=2, title=r"$k^*\sigma$")
    ax_c.tick_params(axis="x", rotation=35)
    polish(ax_c)

    first = runs.groupby("pi").first_mixed_bin.mean()
    ax_d.plot(first.index, first.values, marker="o", color=ORANGE)
    ax_d.set_xlabel(r"$\pi$")
    ax_d.set_ylabel("fraction in first admitted bin")
    ax_d.set_ylim(-0.03, 1.03)
    polish(ax_d)
    for ax, lab in zip(axes.ravel(), "ABCD"):
        panel(ax, lab)
    finish(fig, "FigS3_pi_kappa_replicates")


def figs4_composition_observable_limitation() -> None:
    summ = read("Fig4_composition_summary.csv")
    runs = read("Fig4_composition_runs_symmetric.csv")
    ratios = read("Fig4_composition_symmetry_ratios.csv")
    piv = summ.pivot(index="f_A", columns="kappa", values="C_symmetric_N_mean")
    x, y = piv.columns.to_numpy(float), piv.index.to_numpy(float)
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.9), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()
    im = ax_a.pcolormesh(centers_to_edges(x), centers_to_edges(y), piv.values, shading="flat", cmap="viridis", norm=LogNorm(piv.values.min(), piv.values.max()))
    publication_colorbar(fig, im, ax=ax_a, pad=0.02, label=r"raw $S_{\psi\psi}/2$ peak")
    ax_a.set_xlabel(r"$\kappa$")
    ax_a.set_ylabel(r"$f_A$")

    for kappa, color in ((0.0, BLUE), (0.6, ORANGE), (1.0, PURPLE)):
        sub = summ[np.isclose(summ.kappa, kappa)]
        raw = runs[np.isclose(runs.kappa, kappa)].sort_values(["f_A", "seed"])
        overlay_seed_points(ax_b, raw, "f_A", "C_symmetric_N", color=color, jitter=0.012)
        overlay_seed_points(
            ax_c,
            raw,
            "f_A",
            "C_symmetric_N_comp_normalized",
            color=color,
            jitter=0.012,
        )
        ax_b.errorbar(sub.f_A, sub.C_symmetric_N_mean, yerr=sub.C_symmetric_N_sem, marker="o", ms=3, capsize=2, color=color, label=rf"$\kappa={kappa:g}$", zorder=4)
        ax_c.errorbar(sub.f_A, sub.C_comp_normalized_mean, yerr=sub.C_comp_normalized_sem, marker="o", ms=3, capsize=2, color=color, zorder=4)
    ax_b.set_xlabel(r"$f_A$")
    ax_b.set_ylabel(r"$S_{\psi\psi}^{(N)}(k^*)/2$")
    legend_above(ax_b, ncol=3)
    polish(ax_b)
    ax_c.set_xlabel(r"$f_A$")
    ax_c.set_ylabel(r"$[S_{\psi\psi}^{(N)}/2]/[4f_A(1-f_A)]$")
    polish(ax_c)

    for f_a, color in ((0.2, BLUE), (0.3, GREEN), (0.4, ORANGE)):
        sub = ratios[np.isclose(ratios.f_A, f_a)]
        ax_d.errorbar(sub.kappa, sub.C_ratio_f_over_1minusf, yerr=sub.ratio_propagated_sem, marker="o", ms=3, capsize=2, color=color, label=rf"$f_A={f_a:g}$")
    ax_d.axhline(1, color=BLACK, ls=":", lw=1.0)
    ax_d.set_xlabel(r"$\kappa$")
    ax_d.set_ylabel(r"$S_{\psi\psi}(f)/S_{\psi\psi}(1-f)$")
    legend_above(ax_d, ncol=3)
    polish(ax_d)
    for ax, lab in zip(axes.ravel(), "ABCD"):
        panel(ax, lab)
    finish(fig, "FigS4_composition_observable_limitation")


def figs5_interaction_controls() -> None:
    summ = read("Fig5_interaction_summary.csv")
    runs = read("Fig5_interaction_runs.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.9), constrained_layout=True)
    ax_a, ax_b, ax_c = axes.ravel()
    for seq in ("alternating", "random", "block", "correlated"):
        sub = summ[summ.sequence == seq].sort_values("eps_AB")
        raw = runs[runs.sequence == seq].sort_values(["eps_AB", "seed"])
        overlay_seed_points(
            ax_a,
            raw,
            "eps_AB",
            "C_archive_N",
            color=SEQ_COL[seq],
            jitter=0.08,
            log_x=True,
        )
        ax_a.errorbar(sub.eps_AB, sub.C_archive_N_mean, yerr=sub.C_archive_N_sem, marker="o", ms=2.8, capsize=1.5, color=SEQ_COL[seq], label=seq, zorder=4)
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"$\epsilon_{AB}$")
    ax_a.set_ylabel(r"$\widehat C_N$")
    legend_above(ax_a, ncol=4)
    polish(ax_a)

    for seq in ("alternating", "random", "block", "correlated"):
        sub = summ[summ.sequence == seq].sort_values("eps_AB")
        base = float(sub[np.isclose(sub.eps_AB, 0.1)].C_archive_N_mean.iloc[0])
        ax_b.plot(sub.eps_AB, sub.C_archive_N_mean / base, marker="o", ms=2.8, color=SEQ_COL[seq], label=seq)
    ax_b.set_xscale("log")
    ax_b.set_xlabel(r"$\epsilon_{AB}$")
    ax_b.set_ylabel(r"$\widehat C_N/\widehat C_N(0.1)$")
    polish(ax_b)

    for seq in ("alternating", "random", "block", "correlated"):
        sub = summ[summ.sequence == seq]
        ax_c.plot(sub.eps_AB, sub.k_star_mode, marker="o", ms=2.8, color=SEQ_COL[seq], label=seq)
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"$\epsilon_{AB}$")
    ax_c.set_ylabel(r"modal $k^*\sigma$")
    polish(ax_c)

    for ax, lab in zip(axes.ravel(), "ABC"):
        panel(ax, lab)
    finish(fig, "FigS5_interaction_controls")


def figs6_condition_diagnostics() -> None:
    df = read("FigS6_condition_control_summary.csv").copy()
    runs = read("FigS6_condition_control_runs.csv").copy()
    df["condition"] = df.condition.replace({"dense": "high-density", "soft": "soft tail/T", "short chain": "short-chain"})
    runs["condition"] = runs.condition.replace({"dense": "high-density", "soft": "soft tail/T", "short chain": "short-chain"})
    fig, ax = plt.subplots(figsize=(7.0, 3.0), constrained_layout=True)
    for condition in ("baseline", "high-density", "soft tail/T", "short-chain"):
        sub = df[df.condition == condition]
        raw = runs[runs.condition == condition].sort_values(["kappa", "seed"])
        color = COND_COL[condition]
        overlay_seed_points(ax, raw, "kappa", "T_inst", color=color, jitter=0.018)
        ax.errorbar(sub.kappa, sub.T_inst_mean, yerr=sub.T_inst_sem, marker="o", ms=3, capsize=2, color=color, label=condition, zorder=4)
    ax.axhline(0.7, color=BLACK, ls=":", lw=1.0)
    ax.set_xlabel(r"$\kappa$")
    ax.set_ylabel(r"instantaneous $T^*$")
    legend_above(ax, ncol=4)
    polish(ax)
    finish(fig, "FigS6_condition_diagnostics")


def figs7_rpa_form_factors() -> None:
    form = read("Fig6_finite_chain_form_factor.csv")
    pts = read("Fig6_predictor_points.csv")
    bsens = read("FigS7_predictor_b_sensitivity.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.9), constrained_layout=True)
    ax_a, ax_b, ax_c = axes.ravel()
    for kappa, color in ((0.3, BLUE), (0.6, ORANGE), (1.0, PURPLE)):
        sub = form[(np.isclose(form.pi, 0.99)) & (np.isclose(form.kappa, kappa))]
        ax_a.plot(sub.q_sigma, sub.P_N_q_normalized, color=color, label=rf"$\kappa={kappa:g}$")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"$q\sigma$")
    ax_a.set_ylabel(r"$P_N(q)$")
    legend_above(ax_a, ncol=3)
    polish(ax_a)

    primary = pts[pts.primary_regression_included.astype(bool)]
    sc = ax_b.scatter(primary.inverse_P_N_q0, primary.inverse_C_archive_N, c=primary.pi, cmap="plasma", s=22, edgecolor="white", linewidth=1.0)
    lr = linregress(primary.inverse_P_N_q0, primary.inverse_C_archive_N)
    xx = np.linspace(primary.inverse_P_N_q0.min(), primary.inverse_P_N_q0.max(), 100)
    ax_b.plot(xx, lr.intercept + lr.slope * xx, color=BLACK)
    publication_colorbar(fig, sc, ax=ax_b, pad=0.02, label=r"$\pi$")
    ax_b.set_xlabel(r"$1/P_N(0)$")
    ax_b.set_ylabel(r"$1/\langle\widehat C_N\rangle$")
    polish(ax_b)

    ax_c.plot(bsens.b_over_sigma, bsens.finite_q_r_squared_active_interior, marker="o", color=GREEN, label=r"finite $q=k^*$")
    ax_c.axhline(float(bsens.q0_r_squared_active_interior.iloc[0]), color=GRAY, ls="--", label=r"$q=0$")
    ax_c.axvline(1.0, color=BLACK, ls=":", lw=1.0, label=r"main: $b=\sigma$")
    ax_c.set_xlabel(r"assumed $b/\sigma$")
    ax_c.set_ylim(0.5, 1.0)
    ax_c.set_ylabel(r"$R^2$")
    legend_above(ax_c, ncol=1)
    polish(ax_c)
    for ax, lab in zip(axes.ravel(), "ABC"):
        panel(ax, lab)
    finish(fig, "FigS7_rpa_form_factors")


def main() -> None:
    clean_figure_outputs()
    fig1_sequence_statistics()
    fig2_condition_controls()
    fig3_pi_kappa_response()
    fig4_finite_predictor_rpa()
    fig5_cross_attraction()
    fig6_peak_binning_sensitivity()
    figs1_representative_spectra()
    figs2_sequence_validation()
    figs3_pi_kappa_replicates()
    figs4_composition_observable_limitation()
    figs5_interaction_controls()
    figs6_condition_diagnostics()
    figs7_rpa_form_factors()
    write_submission_aliases()
    sync_journal_destinations()
    print(
        f"Wrote {len(list(PDF.glob('*.pdf')))} source PDFs/PNGs and "
        f"{len(list(SUBMISSION_PDF.glob('*.pdf')))} submission aliases per format; "
        "synchronized journal destinations"
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rpa-only", action="store_true",
                        help="Regenerate only Fig. 4 from existing CSVs; no MD or other figures.")
    parser.add_argument("--source-data", type=Path,
                        help="Directory containing the submitted source_data CSVs.")
    parser.add_argument("--output", type=Path,
                        help="Destination for the RPA-only Fig. 4 PDF and PNG.")
    args = parser.parse_args()
    if args.source_data is not None:
        SD = args.source_data.resolve()
    if args.output is not None and not args.rpa_only:
        parser.error("--output is supported only with --rpa-only")
    if args.rpa_only:
        out = (args.output or HERE / "rpa_correction").resolve()
        out.mkdir(parents=True, exist_ok=True)
        PDF = PNG = out
        fig4_finite_predictor_rpa()
        for suffix in ("pdf", "png"):
            generated = out / f"Fig4_finite_predictor_rpa.{suffix}"
            generated.replace(out / f"Fig4_sequence_response_theory.{suffix}")
        print(f"Wrote corrected Fig. 4 only: {out}")
    else:
        main()
