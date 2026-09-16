"""Compare eps_AB sweep campaigns run at different temperatures.

    python scripts/compare_dsf_campaigns.py --out analysis/dsf_series \
        output/devbox/dsf_T0p7:0.7 output/devbox/dsf_T1p0:1.0 output/devbox/dsf_T1p5:1.5 output/devbox/dsf_T2p0:2.0

Each positional argument is ``<campaign_dir>:<T_star>``. Reads the per-condition, spectra,
F(q*,t) and tau-by-shell tables written by ``melt.epsab_analysis`` and produces:

- ``S_peak_vs_chi_axis.png``: S_psi(q*) against delta_eps / T* for every temperature (the
  Flory-Huggins mapping predicts a single curve when alpha is temperature independent);
- ``F_kmin_vs_lag_by_T.png``: F_psi(q_min, t) at chi = 0 and one mixed-side condition per T;
- ``tau_vs_q_by_T.png``: relaxation time versus wavevector at chi = 0 with q^-2 and q^-4 guides;
- ``tau_kmin_vs_delta_eps.png``: relaxation time of the lowest shell versus delta_eps per T;
- ``dsf_series_summary.csv``: per (T, eps_AB) row of the key numbers.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _f(x: str) -> float:
    return float(x) if x not in ("", "nan", None) else float("nan")


def load(campaign: Path) -> dict:
    a = campaign / "analysis"
    cond = list(csv.DictReader(open(a / "per_condition.csv")))
    fpk = list(csv.DictReader(open(a / "F_peak_by_condition.csv")))
    tau = list(csv.DictReader(open(a / "tau_by_shell_condition.csv")))
    spec = list(csv.DictReader(open(a / "spectra_by_condition.csv")))
    return {"cond": cond, "F": fpk, "tau": tau, "spec": spec}


def main(argv: list[str] | None = None) -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("campaigns", nargs="+", help="<campaign_dir>:<T_star>")
    parser.add_argument("--out", default="analysis/dsf_series")
    parser.add_argument("--eps-like", type=float, default=1.0)
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    series = []
    for item in args.campaigns:
        path, T = item.rsplit(":", 1)
        series.append((float(T), Path(path), load(Path(path))))
    series.sort(key=lambda s: s[0])
    cmap = plt.get_cmap("plasma")
    colors = {T: cmap(0.15 + 0.7 * i / max(len(series) - 1, 1)) for i, (T, _, _) in enumerate(series)}

    # Summary table and the chi-axis collapse.
    rows = []
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for T, path, d in series:
        cond = sorted(d["cond"], key=lambda r: _f(r["delta_eps"]))
        x = np.array([_f(r["delta_eps"]) / T for r in cond])
        y = np.array([_f(r["S_psi_peak_mean"]) for r in cond])
        e = np.array([_f(r["S_psi_peak_sem"]) for r in cond])
        ax.errorbar(x, y, yerr=np.nan_to_num(e), marker="o", ms=4, lw=1.2, capsize=2, color=colors[T], label=f"T*={T:g}")
        for r in cond:
            eps = _f(r["eps_AB"])
            q_min_tau = [_f(t["tau_psi_mean"]) for t in d["tau"] if _f(t["eps_AB"]) == eps]
            rows.append({
                "T_star": T, "eps_AB": eps, "delta_eps": _f(r["delta_eps"]), "delta_eps_over_T": _f(r["delta_eps"]) / T,
                "S_psi_peak_mean": _f(r["S_psi_peak_mean"]), "S_psi_peak_sem": _f(r["S_psi_peak_sem"]),
                "q_peak_mean": _f(r["q_peak_mean"]), "tau_psi_peak_mean": _f(r["tau_psi_peak_mean"]),
                "tau_psi_kmin_mean": q_min_tau[0] if q_min_tau else float("nan"),
                "plateau_psi_peak_mean": _f(r["plateau_psi_peak_mean"]),
                "non_gaussian_ratio_peak_mean": _f(r["non_gaussian_ratio_peak_mean"]),
                "shell_anisotropy_peak_mean": _f(r["shell_anisotropy_peak_mean"]),
                "S_psi_peak_relative_drift_mean": _f(r["S_psi_peak_relative_drift_mean"]),
                "mean_Rg_mean": _f(r["mean_Rg_mean"]),
            })
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta\varepsilon / T^*$  ($\propto\chi$)")
    ax.set_ylabel(r"$S_{\psi\psi}(q^*)$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "S_peak_vs_chi_axis.png", dpi=160)
    fig.savefig(out / "S_peak_vs_chi_axis.pdf")
    plt.close(fig)
    with open(out / "dsf_series_summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # F(q*, t) at chi = 0 and at the first mixed-side step, per temperature.
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
    for T, path, d in series:
        for ax, eps in zip(axes, (1.0, 0.95)):
            rs = [r for r in d["F"] if abs(_f(r["eps_AB"]) - eps) < 1e-9]
            if not rs:
                continue
            t = np.array([_f(r["lag_time"]) for r in rs])
            F = np.array([_f(r["F_psi_mean"]) for r in rs])
            ax.plot(t[1:], F[1:], color=colors[T], lw=1.4, label=f"T*={T:g} (q*={_f(rs[0]['q_peak']):.2f})")
            ax.set_title(rf"$\varepsilon_{{AB}}={eps:g}$  ($\Delta\varepsilon={args.eps_like - eps:g}$)")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel(r"lag $t$ ($\tau$)")
        ax.axhline(1 / np.e, color="0.7", lw=0.8, ls=":")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=7)
    axes[0].set_ylabel(r"$F_{\psi\psi}(q^*,t)$")
    fig.tight_layout()
    fig.savefig(out / "F_kmin_vs_lag_by_T.png", dpi=160)
    fig.savefig(out / "F_kmin_vs_lag_by_T.pdf")
    plt.close(fig)

    # tau(q) at chi = 0 with Rouse/diffusive guides.
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for T, path, d in series:
        rs = [r for r in d["tau"] if abs(_f(r["eps_AB"]) - 1.0) < 1e-9]
        q = np.array([_f(r["shell_q"]) for r in rs])
        tau = np.array([_f(r["tau_psi_mean"]) for r in rs])
        good = np.isfinite(tau) & (tau > 0)
        ax.plot(q[good], tau[good], marker="o", ms=4, lw=1.2, color=colors[T], label=f"T*={T:g}")
    qq = np.array([0.3, 1.2])
    ax.plot(qq, 150.0 * (qq / 0.3) ** -2, "k--", lw=0.8, label=r"$q^{-2}$")
    ax.plot(qq, 150.0 * (qq / 0.3) ** -4, "k:", lw=0.8, label=r"$q^{-4}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$q$ ($\sigma^{-1}$)")
    ax.set_ylabel(r"$\tau_{1/e}(q)$ at $\chi=0$ ($\tau$)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "tau_vs_q_by_T.png", dpi=160)
    fig.savefig(out / "tau_vs_q_by_T.pdf")
    plt.close(fig)

    # Critical slowing down: tau of the lowest few shells versus delta_eps.
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for T, path, d in series:
        qs = sorted({_f(r["shell_q"]) for r in d["tau"]})
        for j, ls in zip((0, 2), ("-", "--")):
            rs = sorted([r for r in d["tau"] if abs(_f(r["shell_q"]) - qs[j]) < 1e-9], key=lambda r: _f(r["delta_eps"]))
            x = np.array([_f(r["delta_eps"]) for r in rs])
            y = np.array([_f(r["tau_psi_mean"]) for r in rs])
            ax.plot(x, y, marker="o", ms=4, lw=1.2, ls=ls, color=colors[T], label=f"T*={T:g}, q={qs[j]:.2f}")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta\varepsilon$")
    ax.set_ylabel(r"$\tau_{1/e}(q)$ ($\tau$); missing = no decay in window")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "tau_kmin_vs_delta_eps.png", dpi=160)
    fig.savefig(out / "tau_kmin_vs_delta_eps.pdf")
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
