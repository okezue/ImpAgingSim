"""Convert a stored ``trajectory.npz`` into ``mode_amplitudes.npz`` for the dynamic structure factor.

Archived campaigns saved per-bead positions at every snapshot but no Fourier amplitudes.
This evaluates the exact box modes from those positions so the archived runs can be
analyzed with ``melt.dynamic_structure`` exactly like newly recorded runs.

    python scripts/modes_from_trajectory.py RUN_DIR [RUN_DIR ...] [--q-max 1.5] [--analyze]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from melt.dynamic_structure import analyze_run  # noqa: E402
from melt.modes import ModeAmplitudeRecorder  # noqa: E402


def convert(run_dir: str, q_max: float, overwrite: bool = False) -> str:
    out = os.path.join(run_dir, "mode_amplitudes.npz")
    if os.path.exists(out) and not overwrite:
        return out
    with np.load(os.path.join(run_dir, "trajectory.npz")) as traj:
        steps = traj["steps"].astype(np.int64)
        positions = traj["positions"]
        types = traj["types"].astype(np.int8)
        box = float(traj["box_size"])
    with open(os.path.join(run_dir, "meta.json")) as handle:
        meta = json.load(handle)
    recorder = ModeAmplitudeRecorder(box, q_max)
    for step, frame in zip(steps, positions):
        recorder.observe(int(step), frame.astype(np.float64), types)
    extra = {
        "mode_interval_steps": int(np.min(np.diff(steps))) if steps.size > 1 else int(steps[0]),
        "production_steps": int(meta["run_params"]["n_steps"]),
        "dt_ps": float(meta["melt_params"]["dt"]),
        "lj_eps_AB": float(meta["melt_params"]["lj_eps_AB"]),
        "seed": int(meta["run_params"]["seed"]),
        "source": "converted from trajectory.npz",
    }
    params = meta.get("sequence_parameters") or {}
    for key in ("kappa", "pi", "f_A"):
        if key in params:
            extra[key] = float(params[key])
    if "T_quench_star" in meta:
        extra["T_quench_star"] = float(meta["T_quench_star"])
    return recorder.save(out, **extra)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dirs", nargs="+")
    parser.add_argument("--q-max", type=float, default=1.5)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--analyze", action="store_true", help="also run melt.dynamic_structure on each run")
    parser.add_argument("--window-fraction", type=float, default=1.0,
                        help="archived runs were annealed beforehand, so the whole production window is stationary")
    args = parser.parse_args(argv)
    for run_dir in args.run_dirs:
        path = convert(run_dir, args.q_max, args.overwrite)
        print(f"{run_dir}: {os.path.basename(path)} written")
        if args.analyze:
            s = analyze_run(run_dir, window_fraction=args.window_fraction)["summary"]
            print(
                f"  q*={s['q_peak']:.3f}  S_psi(q*)={s['S_psi_peak']:.2f}  R4={s['non_gaussian_ratio_peak']:.3f}"
                f"  tau_1/e(q*)={s['tau_psi_peak']:.3g}  plateau={s['plateau_psi_peak']:.3f}"
                f"  drift={s['S_psi_peak_relative_drift']:+.3f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
