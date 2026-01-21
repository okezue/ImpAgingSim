#!/usr/bin/env python3
"""Analyze IMP heteropolymer aging results.

Inputs
------
CSV produced by imp_aging.py.

Outputs
-------
Writes summary CSVs and (optionally) plots:
- summary_per_tw.csv: per condition + waiting time, includes t* and tau
- summary_condition.csv: per condition, includes aging exponent mu and (optional) nu fits

Notes
-----
- Assumes each CSV row is already averaged over trajectories for a single disorder
  realization, as produced by imp_aging.py.
- We compute disorder-averaged curves by averaging over disorder_idx.
- We compute a normalized overlap Q_norm(tw,t) = Q(tw,t)/Q(tw,0) when lag=0 exists.
  This ensures Q_norm(tw,0)=1 and makes tau extraction meaningful.
- chi4(t) is treated as the variance of the overlap across trajectories (already in
  the CSV); we disorder-average it afterward.

The definition chi4 ~ N( <Q^2>-<Q>^2 ) is standard in dynamical heterogeneity work.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class TauResult:
    tau: float
    ok: bool
    reason: str = ""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--csv", type=str, required=True, help="CSV from imp_aging.py")
    p.add_argument("--outdir", type=str, default=".")

    p.add_argument("--target", type=float, default=math.exp(-1), help="Target value for Q_norm to define tau")
    p.add_argument("--mono", action="store_true", help="Monotonize Q_norm by cumulative minimum before crossing")

    # nu fit options
    p.add_argument("--fit_nu", action="store_true", help="Fit short-time scaling exponents nu")
    p.add_argument("--nu_lag_min", type=float, default=3.0)
    p.add_argument("--nu_lag_max", type=float, default=300.0)
    p.add_argument("--nu_tw", type=int, nargs="+", default=[0], help="Which waiting times to fit nu at")

    # plotting
    p.add_argument("--plot", action="store_true", help="Write diagnostic plots")
    p.add_argument("--max_tw_plots", type=int, default=10)

    return p.parse_args()


def disorder_average(df: pd.DataFrame) -> pd.DataFrame:
    gcols = ["ensemble", "epsilon", "tw", "lag"]
    agg_spec = {
        "Q": ("Q_mean", "mean"),
        "chi4": ("chi4", "mean"),
        "D4": ("D4_mean", "mean"),
    }
    if "D2_mean" in df.columns:
        agg_spec["D2"] = ("D2_mean", "mean")
    out = df.groupby(gcols).agg(**agg_spec).reset_index()
    return out


def compute_qnorm(agg: pd.DataFrame) -> pd.DataFrame:
    """Add Q0 and Q_norm columns if lag==0 exists per (ensemble,epsilon,tw)."""
    agg = agg.copy()
    # Map (ensemble,epsilon,tw) -> Q0 at lag=0
    lag0 = agg[agg["lag"] == 0][["ensemble", "epsilon", "tw", "Q"]].rename(columns={"Q": "Q0"})
    if lag0.empty:
        agg["Q0"] = np.nan
        agg["Q_norm"] = np.nan
        return agg

    key = ["ensemble", "epsilon", "tw"]
    agg = agg.merge(lag0, on=key, how="left")
    agg["Q_norm"] = agg["Q"] / agg["Q0"]
    return agg


def first_crossing_tau(lag: np.ndarray, y: np.ndarray, target: float, mono: bool) -> TauResult:
    """Return tau where y crosses down through target.

    Uses log-time interpolation between adjacent points.
    Expects lag sorted ascending.
    """
    if lag.size < 2:
        return TauResult(tau=float("nan"), ok=False, reason="too_few_points")

    # Remove lag==0 for interpolation scale
    mask = lag > 0
    lagp = lag[mask].astype(float)
    yp = y[mask].astype(float)

    if lagp.size < 2:
        return TauResult(tau=float("nan"), ok=False, reason="no_positive_lags")

    if mono:
        yp = np.minimum.accumulate(yp)

    # Need to start above target
    if not (yp[0] >= target):
        return TauResult(tau=float("nan"), ok=False, reason="starts_below_target")

    # If never crosses
    if np.all(yp > target):
        return TauResult(tau=float("nan"), ok=False, reason="no_crossing")

    k = int(np.argmax(yp <= target))
    if k == 0:
        return TauResult(tau=float(lagp[0]), ok=True)

    x0, x1 = math.log(lagp[k - 1]), math.log(lagp[k])
    y0, y1 = float(yp[k - 1]), float(yp[k])
    if y1 == y0:
        return TauResult(tau=float(lagp[k]), ok=True)

    frac = (target - y0) / (y1 - y0)
    frac = min(1.0, max(0.0, frac))
    tau = math.exp(x0 + frac * (x1 - x0))
    return TauResult(tau=tau, ok=True)


def fit_power_law(lag: np.ndarray, y: np.ndarray, lag_min: float, lag_max: float) -> Tuple[float, float, int]:
    """Fit y ~ c * t^nu on a log-log window. Returns (nu, r2, n_pts)."""
    mask = (lag >= lag_min) & (lag <= lag_max) & (lag > 0) & np.isfinite(y) & (y > 0)
    x = np.log(lag[mask].astype(float))
    z = np.log(y[mask].astype(float))
    n = int(x.size)
    if n < 3:
        return float("nan"), float("nan"), n
    nu, logc = np.polyfit(x, z, deg=1)
    # R^2
    zhat = nu * x + logc
    ss_res = float(np.sum((z - zhat) ** 2))
    ss_tot = float(np.sum((z - z.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(nu), float(r2), n


def collapse_score(xs: List[np.ndarray], ys: List[np.ndarray], ngrid: int = 60) -> float:
    """Compute a simple collapse score: average std of y across curves on a common x grid.

    - Interpolates in log x.
    - Uses overlapping x range across curves.
    """
    if len(xs) < 2:
        return float("nan")

    # Determine overlap range
    xmin = max(float(np.min(x[x > 0])) for x in xs)
    xmax = min(float(np.max(x)) for x in xs)
    if not (xmax > xmin > 0):
        return float("nan")

    grid = np.logspace(math.log10(xmin), math.log10(xmax), ngrid)

    ygrid = []
    for x, y in zip(xs, ys):
        mask = (x > 0) & np.isfinite(y)
        x1 = x[mask].astype(float)
        y1 = y[mask].astype(float)
        if x1.size < 2:
            return float("nan")
        # Interpolate in log space
        lx = np.log(x1)
        ly = y1
        lgrid = np.log(grid)
        yi = np.interp(lgrid, lx, ly)
        ygrid.append(yi)

    Y = np.stack(ygrid, axis=0)
    return float(np.mean(np.std(Y, axis=0)))


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    if df.empty:
        raise SystemExit("Empty CSV")

    agg = disorder_average(df)
    agg = compute_qnorm(agg)

    # Per-(ensemble,epsilon,tw) summaries
    per_tw_rows: List[Dict[str, float]] = []
    cond_rows: List[Dict[str, float]] = []

    for (ens, eps), subc in agg.groupby(["ensemble", "epsilon"]):
        subc = subc.sort_values(["tw", "lag"]).copy()

        # Compute per-tw t* and tau
        taus: List[Tuple[int, float]] = []
        for tw, sub in subc.groupby("tw"):
            sub = sub.sort_values("lag")
            lag = sub["lag"].to_numpy(dtype=float)
            chi4 = sub["chi4"].to_numpy(dtype=float)

            # t* excludes lag=0 if present
            mask_pos = lag > 0
            if np.any(mask_pos):
                idx = int(np.argmax(chi4[mask_pos]))
                lpos = lag[mask_pos]
                cpos = chi4[mask_pos]
                t_star = float(lpos[idx])
                chi4_star = float(cpos[idx])
            else:
                t_star = float("nan")
                chi4_star = float("nan")

            Q0 = float(sub[sub["lag"] == 0]["Q0"].iloc[0]) if np.any(sub["lag"] == 0) and "Q0" in sub.columns else float("nan")
            Qn = sub["Q_norm"].to_numpy(dtype=float) if "Q_norm" in sub.columns else np.full_like(lag, np.nan)
            tau_res = first_crossing_tau(lag, Qn, target=float(args.target), mono=bool(args.mono))
            tau = float(tau_res.tau)

            if tw > 0 and math.isfinite(tau):
                taus.append((int(tw), tau))

            row = {
                "ensemble": ens,
                "epsilon": float(eps),
                "tw": float(tw),
                "t_star": t_star,
                "chi4_star": chi4_star,
                "Q0": Q0,
                "tau": tau,
            }
            per_tw_rows.append(row)

        # Fit aging exponent mu from tau(tw) ~ tw^mu
        mu = float("nan")
        mu_se = float("nan")
        nfit = 0
        if len(taus) >= 2:
            tws = np.array([t for t, _ in taus], dtype=float)
            tauv = np.array([x for _, x in taus], dtype=float)
            x = np.log(tws)
            y = np.log(tauv)
            nfit = int(x.size)
            # linear regression
            A = np.vstack([x, np.ones_like(x)]).T
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            mu = float(coef[0])
            # standard error
            yhat = A @ coef
            resid = y - yhat
            s2 = float(np.sum(resid**2) / max(1, nfit - 2))
            xtx_inv = np.linalg.inv(A.T @ A)
            mu_se = float(math.sqrt(s2 * xtx_inv[0, 0])) if nfit > 2 else float("nan")

        cond_row: Dict[str, float] = {
            "ensemble": ens,
            "epsilon": float(eps),
            "mu": mu,
            "mu_se": mu_se,
            "mu_n": float(nfit),
        }

        # Optional nu fits (short-time) on disorder-averaged curves
        if args.fit_nu:
            for tw_fit in args.nu_tw:
                sub = subc[subc["tw"] == tw_fit].sort_values("lag")
                if sub.empty:
                    continue
                lag = sub["lag"].to_numpy(dtype=float)
                if "D2" in sub.columns:
                    nu2, r2_2, n2 = fit_power_law(lag, sub["D2"].to_numpy(dtype=float), args.nu_lag_min, args.nu_lag_max)
                    cond_row[f"nu_D2_tw{tw_fit}"] = nu2
                    cond_row[f"nu_D2_r2_tw{tw_fit}"] = r2_2
                    cond_row[f"nu_D2_n_tw{tw_fit}"] = float(n2)
                nu4, r2_4, n4 = fit_power_law(lag, sub["D4"].to_numpy(dtype=float), args.nu_lag_min, args.nu_lag_max)
                cond_row[f"nu_D4_tw{tw_fit}"] = nu4
                cond_row[f"nu_D4_r2_tw{tw_fit}"] = r2_4
                cond_row[f"nu_D4_n_tw{tw_fit}"] = float(n4)

        # Collapse scores: Q_norm vs t/tau and chi4 vs t/t*
        # Use only tw with finite tau and finite t_star
        xs_q: List[np.ndarray] = []
        ys_q: List[np.ndarray] = []
        xs_c: List[np.ndarray] = []
        ys_c: List[np.ndarray] = []

        for tw, sub in subc.groupby("tw"):
            sub = sub.sort_values("lag")
            lag = sub["lag"].to_numpy(dtype=float)
            chi4 = sub["chi4"].to_numpy(dtype=float)
            Qn = sub["Q_norm"].to_numpy(dtype=float) if "Q_norm" in sub.columns else None

            # tau
            tau_vals = [r["tau"] for r in per_tw_rows if r["ensemble"] == ens and r["epsilon"] == float(eps) and r["tw"] == float(tw)]
            tau = float(tau_vals[0]) if tau_vals else float("nan")

            if Qn is not None and math.isfinite(tau) and tau > 0:
                xs_q.append(lag / tau)
                ys_q.append(Qn)

            # t_star and chi4_star
            t_star_vals = [r["t_star"] for r in per_tw_rows if r["ensemble"] == ens and r["epsilon"] == float(eps) and r["tw"] == float(tw)]
            chi4_star_vals = [r["chi4_star"] for r in per_tw_rows if r["ensemble"] == ens and r["epsilon"] == float(eps) and r["tw"] == float(tw)]
            t_star = float(t_star_vals[0]) if t_star_vals else float("nan")
            chi4_star = float(chi4_star_vals[0]) if chi4_star_vals else float("nan")
            if math.isfinite(t_star) and t_star > 0 and math.isfinite(chi4_star) and chi4_star > 0:
                xs_c.append(lag / t_star)
                ys_c.append(chi4 / chi4_star)

        cond_row["collapse_Qnorm"] = collapse_score(xs_q, ys_q) if len(xs_q) >= 2 else float("nan")
        cond_row["collapse_chi4"] = collapse_score(xs_c, ys_c) if len(xs_c) >= 2 else float("nan")

        cond_rows.append(cond_row)

        # Optional plots
        if args.plot:
            import matplotlib.pyplot as plt

            # Limit waiting times plotted
            tws_sorted = sorted(subc["tw"].unique())[: int(args.max_tw_plots)]

            # Q_norm vs lag and collapse
            if "Q_norm" in subc.columns and np.any(subc["lag"] == 0):
                plt.figure()
                for tw in tws_sorted:
                    s = subc[subc["tw"] == tw].sort_values("lag")
                    plt.plot(s["lag"], s["Q_norm"], marker="o", linestyle="-")
                plt.xscale("log")
                plt.ylim(0, 1.05)
                plt.xlabel("lag (MC sweeps)")
                plt.ylabel("Q_norm")
                plt.title(f"Q_norm vs lag | ens={ens} eps={eps}")
                plt.tight_layout()
                plt.savefig(os.path.join(args.outdir, f"Qnorm_lag_ens-{ens}_eps-{eps:g}.png"), dpi=200)
                plt.close()

                # Collapsed
                plt.figure()
                for tw in tws_sorted:
                    s = subc[subc["tw"] == tw].sort_values("lag")
                    lag = s["lag"].to_numpy(dtype=float)
                    tau_vals = [r["tau"] for r in per_tw_rows if r["ensemble"] == ens and r["epsilon"] == float(eps) and r["tw"] == float(tw)]
                    tau = float(tau_vals[0]) if tau_vals else float("nan")
                    if math.isfinite(tau) and tau > 0:
                        plt.plot(lag / tau, s["Q_norm"], marker="o", linestyle="-")
                plt.xscale("log")
                plt.ylim(0, 1.05)
                plt.xlabel("lag / tau(tw)")
                plt.ylabel("Q_norm")
                plt.title(f"Q_norm collapse | ens={ens} eps={eps}")
                plt.tight_layout()
                plt.savefig(os.path.join(args.outdir, f"Qnorm_collapse_ens-{ens}_eps-{eps:g}.png"), dpi=200)
                plt.close()

            # chi4 and collapse
            plt.figure()
            for tw in tws_sorted:
                s = subc[subc["tw"] == tw].sort_values("lag")
                plt.plot(s["lag"], s["chi4"], marker="o", linestyle="-")
            plt.xscale("log")
            plt.xlabel("lag (MC sweeps)")
            plt.ylabel("chi4")
            plt.title(f"chi4 vs lag | ens={ens} eps={eps}")
            plt.tight_layout()
            plt.savefig(os.path.join(args.outdir, f"chi4_lag_ens-{ens}_eps-{eps:g}.png"), dpi=200)
            plt.close()

            plt.figure()
            for tw in tws_sorted:
                s = subc[subc["tw"] == tw].sort_values("lag")
                lag = s["lag"].to_numpy(dtype=float)
                chi4 = s["chi4"].to_numpy(dtype=float)
                # t* and chi4*
                mask_pos = lag > 0
                if np.any(mask_pos):
                    idx = int(np.argmax(chi4[mask_pos]))
                    t_star = float(lag[mask_pos][idx])
                    chi4_star = float(chi4[mask_pos][idx])
                    if t_star > 0 and chi4_star > 0:
                        plt.plot(lag / t_star, chi4 / chi4_star, marker="o", linestyle="-")
            plt.xscale("log")
            plt.xlabel("lag / t*(tw)")
            plt.ylabel("chi4 / chi4*")
            plt.title(f"chi4 collapse | ens={ens} eps={eps}")
            plt.tight_layout()
            plt.savefig(os.path.join(args.outdir, f"chi4_collapse_ens-{ens}_eps-{eps:g}.png"), dpi=200)
            plt.close()

    per_tw = pd.DataFrame(per_tw_rows)
    cond = pd.DataFrame(cond_rows)

    per_tw_path = os.path.join(args.outdir, "summary_per_tw.csv")
    cond_path = os.path.join(args.outdir, "summary_condition.csv")

    per_tw.to_csv(per_tw_path, index=False)
    cond.to_csv(cond_path, index=False)

    print(f"Wrote: {per_tw_path}")
    print(f"Wrote: {cond_path}")


if __name__ == "__main__":
    main()
