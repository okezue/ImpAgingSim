"""Phase boundary of the correlated random copolymer melt: crossover eps_AB versus kappa.

    python scripts/compare_kappa_campaigns.py --out analysis/kappa_boundary --T-star 1.5 --alpha 2.0 \
        output/devbox/kappa0p0_T1p5:0.0 output/devbox/kappa0p25_T1p5:0.25 ...

Each positional argument is ``<campaign_dir>:<kappa>``.  For every kappa the script reads the
per-condition table written by ``melt.epsab_analysis`` and locates the mixed -> demixed
crossover with three estimators (half-rise of S_psi(q*), the 1.5 crossing of the per-mode
non-Gaussianity ratio, and the maximum of the peak-intensity time variance), then compares
with the mean-field spinodal of ``melt.rpa`` for the same sequence statistics using the
contact factor ``alpha`` fitted in Stage 2.  Outputs a table, the boundary figure, the
S_psi(q*) curves and the shell anisotropy (lamellar indicator) per kappa.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from melt.rpa import eps_AB_from_chi, rpa_spinodal_chi, segment_length_from_rg  # noqa: E402


def _f(x) -> float:
    return float(x) if x not in ("", "nan", None) else float("nan")


def half_rise(x: np.ndarray, y: np.ndarray, plateau_from: float = 0.3) -> float:
    good = np.isfinite(y)
    x, y = x[good], y[good]
    if x.size < 3:
        return float("nan")
    base = y[0]
    tail = y[x >= plateau_from]
    plateau = tail.mean() if tail.size else y[-3:].mean()
    if plateau <= 1.5 * base:
        return float("nan")
    level = 0.5 * (base + plateau)
    idx = np.where(y >= level)[0]
    if idx.size == 0 or idx[0] == 0:
        return float("nan")
    i = idx[0]
    return float(x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]))


def crossing(x: np.ndarray, y: np.ndarray, level: float) -> float:
    good = np.isfinite(y)
    x, y = x[good], y[good]
    for i in range(x.size - 1):
        if (y[i] - level) * (y[i + 1] - level) <= 0 and y[i] != y[i + 1]:
            return float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return float("nan")


def main(argv: list[str] | None = None) -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("campaigns", nargs="+", help="<campaign_dir>:<kappa>")
    parser.add_argument("--out", default="analysis/kappa_boundary")
    parser.add_argument("--T-star", dest="T_star", type=float, default=1.5)
    parser.add_argument("--alpha", type=float, default=2.0, help="contact factor chi = alpha*delta_eps/T*")
    parser.add_argument("--pi", type=float, default=0.99)
    parser.add_argument("--chain-length", type=int, default=40)
    parser.add_argument("--eps-like", type=float, default=1.0)
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    series = []
    for item in args.campaigns:
        path, kappa = item.rsplit(":", 1)
        cond = list(csv.DictReader(open(Path(path) / "analysis" / "per_condition.csv")))
        spec = list(csv.DictReader(open(Path(path) / "analysis" / "spectra_by_condition.csv")))
        series.append((float(kappa), sorted(cond, key=lambda r: _f(r["delta_eps"])), spec))
    series.sort(key=lambda s: s[0])
    cmap = plt.get_cmap("viridis")
    colors = {k: cmap(0.1 + 0.8 * i / max(len(series) - 1, 1)) for i, (k, _, _) in enumerate(series)}

    rows = []
    fig_s, ax_s = plt.subplots(figsize=(6.0, 4.2))
    fig_a, ax_a = plt.subplots(figsize=(6.0, 4.2))
    fig_r, ax_r = plt.subplots(figsize=(6.0, 4.2))
    for kappa, cond, spec in series:
        x = np.array([_f(r["delta_eps"]) for r in cond])
        S = np.array([_f(r["S_psi_peak_mean"]) for r in cond])
        Se = np.array([_f(r["S_psi_peak_sem"]) for r in cond])
        R4 = np.array([_f(r["non_gaussian_ratio_peak_mean"]) for r in cond])
        aniso = np.array([_f(r["shell_anisotropy_peak_mean"]) for r in cond])
        relvar = np.array([_f(r["peak_intensity_relative_variance_mean"]) for r in cond])
        rg = _f(cond[0]["mean_Rg_mean"])
        b = segment_length_from_rg(rg, args.chain_length) if np.isfinite(rg) and rg > 0 else 1.0
        L = _f(cond[0]["box_size"])
        chi_s, _ = rpa_spinodal_chi(args.chain_length, kappa, args.pi, 0.5, b, q=np.asarray([2 * np.pi / L]))
        eps_s = float(eps_AB_from_chi(chi_s, args.alpha, args.T_star, args.eps_like))
        d_half = half_rise(x, S)
        d_r4 = crossing(x, R4, 1.5)
        d_var = float(x[np.nanargmax(relvar)]) if np.any(np.isfinite(relvar)) else float("nan")
        rows.append({
            "kappa": kappa, "T_star": args.T_star, "Rg_chi0": rg, "segment_length_b": b,
            "S_peak_chi0": S[0], "S_peak_max": np.nanmax(S),
            "delta_eps_half_rise": d_half, "eps_AB_half_rise": args.eps_like - d_half if np.isfinite(d_half) else float("nan"),
            "delta_eps_R4_1p5": d_r4, "eps_AB_R4_1p5": args.eps_like - d_r4 if np.isfinite(d_r4) else float("nan"),
            "delta_eps_max_intensity_variance": d_var,
            "chi_spinodal_rpa": float(chi_s), "delta_eps_spinodal_rpa": args.eps_like - eps_s, "eps_AB_spinodal_rpa": eps_s,
            "max_shell_anisotropy": np.nanmax(aniso), "delta_eps_at_max_anisotropy": float(x[np.nanargmax(aniso)]),
        })
        label = f"kappa={kappa:g}"
        ax_s.errorbar(x, S, yerr=np.nan_to_num(Se), marker="o", ms=4, lw=1.2, capsize=2, color=colors[kappa], label=label)
        ax_a.plot(x, aniso, marker="o", ms=4, lw=1.2, color=colors[kappa], label=label)
        ax_r.plot(x, R4, marker="o", ms=4, lw=1.2, color=colors[kappa], label=label)

    ax_s.set_yscale("log")
    ax_s.set_xlabel(r"$\Delta\varepsilon$")
    ax_s.set_ylabel(rf"$S_{{\psi\psi}}(q^*)$ at $T^*={args.T_star:g}$")
    ax_s.legend(fontsize=8)
    fig_s.tight_layout(); fig_s.savefig(out / "S_peak_vs_delta_eps_by_kappa.png", dpi=160); fig_s.savefig(out / "S_peak_vs_delta_eps_by_kappa.pdf"); plt.close(fig_s)
    ax_a.axhline(3.0, color="0.6", lw=0.8, ls=":"); ax_a.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax_a.set_xlabel(r"$\Delta\varepsilon$"); ax_a.set_ylabel("shell anisotropy at $q^*$ (3 = single lamellar direction)")
    ax_a.legend(fontsize=8)
    fig_a.tight_layout(); fig_a.savefig(out / "anisotropy_vs_delta_eps_by_kappa.png", dpi=160); fig_a.savefig(out / "anisotropy_vs_delta_eps_by_kappa.pdf"); plt.close(fig_a)
    ax_r.axhline(2.0, color="0.6", lw=0.8, ls=":"); ax_r.axhline(1.0, color="0.6", lw=0.8, ls=":"); ax_r.axhline(1.5, color="C3", lw=0.8, ls="--")
    ax_r.set_xlabel(r"$\Delta\varepsilon$"); ax_r.set_ylabel(r"$\langle|\rho_\psi|^4\rangle_t/\langle|\rho_\psi|^2\rangle_t^2$ at $q^*$")
    ax_r.legend(fontsize=8)
    fig_r.tight_layout(); fig_r.savefig(out / "non_gaussian_ratio_vs_delta_eps_by_kappa.png", dpi=160); fig_r.savefig(out / "non_gaussian_ratio_vs_delta_eps_by_kappa.pdf"); plt.close(fig_r)

    # Phase boundary.
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    k = np.array([r["kappa"] for r in rows])
    for key, marker, label in (
        ("delta_eps_half_rise", "o", r"half-rise of $S_{\psi\psi}(q^*)$"),
        ("delta_eps_R4_1p5", "s", "non-Gaussianity ratio = 1.5"),
        ("delta_eps_max_intensity_variance", "^", "max intensity variance"),
    ):
        y = np.array([r[key] for r in rows])
        ax.plot(k, y, marker=marker, ms=6, lw=1.0, label=label)
    ax.plot(k, [r["delta_eps_spinodal_rpa"] for r in rows], "k--", lw=1.2, label=rf"RPA spinodal ($\alpha={args.alpha:g}$)")
    ax.set_xlabel(r"sequence correlation $\kappa$ ($\pi=0.99$)")
    ax.set_ylabel(rf"crossover $\Delta\varepsilon_c$ at $T^*={args.T_star:g}$  (mixed below, demixed above)")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "phase_boundary_kappa.png", dpi=160); fig.savefig(out / "phase_boundary_kappa.pdf"); plt.close(fig)

    with open(out / "phase_boundary_kappa.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow({key: ("" if isinstance(v, float) and not np.isfinite(v) else v) for key, v in r.items()})
    for r in rows:
        print(f"kappa={r['kappa']:g}: half-rise delta_eps={r['delta_eps_half_rise']:.3f}  R4=1.5 at {r['delta_eps_R4_1p5']:.3f}  "
              f"RPA spinodal {r['delta_eps_spinodal_rpa']:.3f}  max anisotropy {r['max_shell_anisotropy']:.2f}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
