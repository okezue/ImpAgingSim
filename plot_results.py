#!/usr/bin/env python3
"""Quick plotting helper for imp_aging.py results.

Reads the CSV produced by imp_aging.py and produces disorder-averaged
curves of Q_mean and chi4 versus lag for each waiting time t_w.

Example:
  python3 plot_results.py --csv results.csv --ensemble correlated --epsilon 6

This script is optional; it is just a convenience for exploratory work.
"""

from __future__ import annotations

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--csv", type=str, required=True)
    p.add_argument("--ensemble", type=str, default="correlated")
    p.add_argument("--epsilon", type=float, default=6.0)
    p.add_argument("--outdir", type=str, default=".")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.csv)
    df = df[(df["ensemble"] == args.ensemble) & (df["epsilon"] == args.epsilon)].copy()
    if df.empty:
        raise SystemExit("No rows matched the requested filters.")

    # disorder average (each row is already averaged over trajectories)
    gcols = ["ensemble", "epsilon", "tw", "lag"]
    agg = (
        df.groupby(gcols)
        .agg(Q_mean=("Q_mean", "mean"), chi4=("chi4", "mean"), D4_mean=("D4_mean", "mean"))
        .reset_index()
    )

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    for tw, sub in agg.groupby("tw"):
        sub = sub.sort_values("lag")

        # Q
        plt.figure()
        plt.plot(sub["lag"], sub["Q_mean"], marker="o", linestyle="-")
        plt.xscale("log")
        plt.xlabel("lag (MC sweeps)")
        plt.ylabel("Q_mean")
        plt.title(f"Q vs lag | ensemble={args.ensemble} eps={args.epsilon} tw={tw}")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"Q_ens-{args.ensemble}_eps-{args.epsilon:g}_tw-{tw}.png"), dpi=200)
        plt.close()

        # chi4
        plt.figure()
        plt.plot(sub["lag"], sub["chi4"], marker="o", linestyle="-")
        plt.xscale("log")
        plt.xlabel("lag (MC sweeps)")
        plt.ylabel("chi4")
        plt.title(f"chi4 vs lag | ensemble={args.ensemble} eps={args.epsilon} tw={tw}")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"chi4_ens-{args.ensemble}_eps-{args.epsilon:g}_tw-{tw}.png"), dpi=200)
        plt.close()

        # D4
        plt.figure()
        plt.plot(sub["lag"], sub["D4_mean"], marker="o", linestyle="-")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("lag (MC sweeps)")
        plt.ylabel("D4_mean")
        plt.title(f"D4 vs lag | ensemble={args.ensemble} eps={args.epsilon} tw={tw}")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"D4_ens-{args.ensemble}_eps-{args.epsilon:g}_tw-{tw}.png"), dpi=200)
        plt.close()

    print(f"Saved plots to: {outdir}")


if __name__ == "__main__":
    main()
