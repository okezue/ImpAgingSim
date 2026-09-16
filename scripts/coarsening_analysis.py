"""Time-resolved coarsening of the microphase pattern from recorded box modes.

    python scripts/coarsening_analysis.py --out analysis/coarsening CAMPAIGN_DIR [--n-bins 40]

For every run in the campaign the composition spectrum S_psi(q, t) is formed frame by
frame from ``mode_amplitudes.npz``, then averaged in logarithmically spaced time bins.
Reported per (eps_AB, time bin), mean and SEM over seeds: the peak height S_psi(q*, t),
the peak wavevector q*(t) (intensity-weighted over the three strongest shells so it moves
continuously between the discrete box shells), the lowest-shell intensity S_psi(q_min, t)
and the coarse-grained variance.  A power-law fit of S_psi(q*, t) over the last decade
gives the growth exponent.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from melt.dynamic_structure import composition_amplitudes, equal_time_spectrum, shell_average_series  # noqa: E402
from melt.modes import load_mode_amplitudes  # noqa: E402


def run_series(run_dir: Path) -> dict[str, np.ndarray]:
    data = load_mode_amplitudes(str(run_dir / "mode_amplitudes.npz"))
    steps = data["steps"].astype(np.int64)
    dt = float(data["dt_ps"]) if "dt_ps" in data else 1.0
    n_total = int(data["n_total"])
    shell_index = data["shell_index"].astype(np.int64)
    shell_q = data["shell_q"].astype(np.float64)
    psi = composition_amplitudes(data)
    S_modes = equal_time_spectrum(psi, n_total)
    S_shell = shell_average_series(S_modes, shell_index, shell_q.size)  # (T, S)
    peak = np.argmax(S_shell, axis=1)
    # Intensity-weighted q over the strongest three shells per frame.
    order = np.argsort(S_shell, axis=1)[:, -3:]
    w = np.take_along_axis(S_shell, order, axis=1)
    q_w = (w * shell_q[order]).sum(axis=1) / w.sum(axis=1)
    return {
        "time": steps * dt,
        "S_peak": S_shell[np.arange(S_shell.shape[0]), peak],
        "q_peak_shell": shell_q[peak],
        "q_peak_weighted": q_w,
        "S_kmin": S_shell[:, 0],
        "eps_AB": float(data["lj_eps_AB"]),
        "seed": int(data["seed"]),
    }


def log_bins(t: np.ndarray, n_bins: int) -> np.ndarray:
    t_pos = t[t > 0]
    return np.geomspace(t_pos.min(), t_pos.max() * 1.0001, n_bins + 1)


def main(argv: list[str] | None = None) -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("campaign")
    parser.add_argument("--out", default="analysis/coarsening")
    parser.add_argument("--n-bins", type=int, default=40)
    args = parser.parse_args(argv)
    campaign = Path(args.campaign)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((campaign / "manifest.json").read_text())

    runs = []
    for spec in manifest["runs"]:
        run_dir = campaign / "runs" / spec["run_id"]
        if (run_dir / "mode_amplitudes.npz").is_file():
            runs.append(run_series(run_dir))
    if not runs:
        raise SystemExit("no runs with mode amplitudes")
    edges = log_bins(runs[0]["time"], args.n_bins)
    centers = np.sqrt(edges[:-1] * edges[1:])
    eps_values = sorted({r["eps_AB"] for r in runs}, reverse=True)

    rows = []
    fits = {}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))
    cmap = plt.get_cmap("viridis")
    for i, eps in enumerate(eps_values):
        group = [r for r in runs if r["eps_AB"] == eps]
        binned = {k: [] for k in ("S_peak", "q_peak_weighted", "S_kmin")}
        for r in group:
            idx = np.digitize(r["time"], edges) - 1
            for k in binned:
                vals = np.full(centers.size, np.nan)
                for b in range(centers.size):
                    sel = idx == b
                    if np.any(sel):
                        vals[b] = np.nanmean(r[k][sel])
                binned[k].append(vals)
        stats = {}
        for k, arr in binned.items():
            a = np.array(arr)
            stats[k] = (np.nanmean(a, axis=0), np.nanstd(a, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(a), axis=0)))
        for b in range(centers.size):
            rows.append({"eps_AB": eps, "delta_eps": 1.0 - eps, "time": centers[b],
                         "S_peak_mean": stats["S_peak"][0][b], "S_peak_sem": stats["S_peak"][1][b],
                         "q_peak_mean": stats["q_peak_weighted"][0][b], "q_peak_sem": stats["q_peak_weighted"][1][b],
                         "S_kmin_mean": stats["S_kmin"][0][b], "S_kmin_sem": stats["S_kmin"][1][b], "n_seeds": len(group)})
        # Growth exponent over the last decade.
        t = centers
        y = stats["S_peak"][0]
        sel = np.isfinite(y) & (t >= t[np.isfinite(y)].max() / 10.0)
        if np.count_nonzero(sel) >= 4:
            slope, intercept = np.polyfit(np.log(t[sel]), np.log(y[sel]), 1)
            yq = stats["q_peak_weighted"][0]
            selq = np.isfinite(yq) & sel
            slope_q = np.polyfit(np.log(t[selq]), np.log(yq[selq]), 1)[0] if np.count_nonzero(selq) >= 4 else float("nan")
            fits[f"{eps:g}"] = {"S_peak_growth_exponent_last_decade": float(slope), "q_peak_exponent_last_decade": float(slope_q),
                                "t_range": [float(t[sel].min()), float(t[sel].max())]}
        color = cmap(0.15 + 0.7 * i / max(len(eps_values) - 1, 1))
        label = rf"$\varepsilon_{{AB}}={eps:g}$"
        axes[0].errorbar(centers, stats["S_peak"][0], yerr=np.nan_to_num(stats["S_peak"][1]), color=color, lw=1.2, label=label)
        axes[1].errorbar(centers, stats["q_peak_weighted"][0], yerr=np.nan_to_num(stats["q_peak_weighted"][1]), color=color, lw=1.2, label=label)
        axes[2].errorbar(centers, stats["S_kmin"][0], yerr=np.nan_to_num(stats["S_kmin"][1]), color=color, lw=1.2, label=label)
    for ax, ylabel in zip(axes, (r"$S_{\psi\psi}(q^*,t)$", r"$q^*(t)$ (intensity-weighted, $\sigma^{-1}$)", r"$S_{\psi\psi}(q_{\min},t)$")):
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"time after quench $t$ ($\tau$)"); ax.set_ylabel(ylabel)
    axes[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "coarsening.png", dpi=160); fig.savefig(out / "coarsening.pdf"); plt.close(fig)

    with open(out / "coarsening_binned.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: ("" if isinstance(v, float) and not np.isfinite(v) else v) for k, v in r.items()})
    with open(out / "coarsening_fits.json", "w") as handle:
        json.dump(fits, handle, indent=2)
    for eps, f in fits.items():
        print(f"eps_AB={eps}: S_peak ~ t^{f['S_peak_growth_exponent_last_decade']:+.3f}, q* ~ t^{f['q_peak_exponent_last_decade']:+.3f} over t in {f['t_range'][0]:.0f}-{f['t_range'][1]:.0f} tau")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
