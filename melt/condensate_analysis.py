"""Aggregate a condensation campaign (protein or chromatin mode).

Per run: trailing-window means of the condensation observables, the B-density spectrum
and its two-time correlation ``F_BB(q*, tau)`` (blob lifetime), and for dynamic marks the
mark autocorrelation (memory time), the persistence fraction and the B-fraction
trajectory.  Per condition: mean and SEM over seeds, phase-diagram heatmaps over the two
varied axes, and time-series figures.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .campaign import campaign_commit, completion_valid, validate_manifest
from .dynamic_structure import (
    block_sem,
    equal_time_spectrum,
    relaxation_times,
    shell_average_series,
    stationary_window_mask,
    time_correlation_by_shell,
)
from .marks import TYPE_B, mark_autocorrelation, persistence_fraction
from .modes import load_mode_amplitudes

DESIGN_KEYS = ("schema_version", "study", "campaign_id", "design", "simulation_settings", "runs")
COND_MEANS = ("f_B", "largest_B_cluster_fraction", "n_B_clusters", "dense_B_fraction", "mean_BB_coordination",
              "mean_A_coordination", "mean_B_coordination", "B_cluster_size_second_moment", "mean_coordination")
RUN_KEYS = ("eps_BB", "kappa", "f_A", "k_off", "k_on", "k_fb", "seed")
OBSERVABLES = COND_MEANS + ("mean_Rg", "S_BB_peak", "q_BB_peak", "S_BB_kmin", "tau_BB_peak", "plateau_BB_peak",
                           "S_psi_peak", "f_B_final", "f_B_drift", "mark_memory_time", "mark_persistence",
                           "largest_cluster_time_std", "f_B_time_std")


def _f(x: Any) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _mean_sem(values: list[float]) -> tuple[float, float]:
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan")
    return float(x.mean()), float(x.std(ddof=1) / np.sqrt(x.size)) if x.size > 1 else float("nan")


def condensation_series(run_dir: Path) -> dict[str, np.ndarray]:
    rows = list(csv.DictReader(open(run_dir / "condensation.csv")))
    return {k: np.asarray([_f(r[k]) for r in rows]) for k in rows[0]}


def analyze_run(run_dir: Path, window_fraction: float, dt: float, production_steps: int) -> dict[str, Any]:
    out: dict[str, Any] = {}
    cond = condensation_series(run_dir)
    steps = cond["step"].astype(np.int64)
    window = steps > (1.0 - window_fraction) * production_steps
    for key in COND_MEANS:
        out[key] = float(np.nanmean(cond[key][window])) if key in cond else float("nan")
    out["largest_cluster_time_std"] = float(np.nanstd(cond["largest_B_cluster_fraction"][window]))
    out["f_B_time_std"] = float(np.nanstd(cond["f_B"][window]))
    out["f_B_final"] = float(cond["f_B"][-1])
    half = np.count_nonzero(window) // 2
    fb_w = cond["f_B"][window]
    out["f_B_drift"] = float(fb_w[half:].mean() - fb_w[:half].mean()) if half > 0 else float("nan")

    snaps = list(csv.DictReader(open(run_dir / "snapshots.csv")))
    rg = [_f(r["mean_Rg"]) for r in snaps if int(r["step"]) > (1.0 - window_fraction) * production_steps]
    out["mean_Rg"] = float(np.mean(rg)) if rg else float("nan")

    data = load_mode_amplitudes(str(run_dir / "mode_amplitudes.npz"))
    msteps = data["steps"].astype(np.int64)
    mask = stationary_window_mask(msteps, production_steps, window_fraction)
    rho_B = np.asarray(data["rho_B"], dtype=np.complex128)[mask]
    psi = np.asarray(data["rho_A"], dtype=np.complex128)[mask] - rho_B
    n_total = int(data["n_total"])
    shell_index = data["shell_index"].astype(np.int64)
    shell_q = data["shell_q"].astype(np.float64)
    S_BB_t = shell_average_series(equal_time_spectrum(rho_B, n_total), shell_index, shell_q.size)
    S_BB = S_BB_t.mean(axis=0)
    peak = int(np.argmax(S_BB))
    out["S_BB_peak"] = float(S_BB[peak])
    out["q_BB_peak"] = float(shell_q[peak])
    out["S_BB_kmin"] = float(S_BB[0])
    S_psi = shell_average_series(equal_time_spectrum(psi, n_total), shell_index, shell_q.size).mean(axis=0)
    out["S_psi_peak"] = float(S_psi.max())
    T = rho_B.shape[0]
    max_lag = max(1, int(0.5 * (T - 1)))
    frame = int(np.min(np.diff(msteps[mask]))) if T > 1 else 1
    S_qt = time_correlation_by_shell(rho_B, shell_index, shell_q.size, n_total, max_lag=max_lag)
    with np.errstate(invalid="ignore", divide="ignore"):
        F = S_qt / S_qt[:, :1]
    lag_time = np.arange(max_lag + 1) * frame * dt
    out["tau_BB_peak"] = float(relaxation_times(lag_time, F)[peak])
    out["plateau_BB_peak"] = float(F[peak, max(1, max_lag // 2):].mean())
    arrays = {"shell_q": shell_q, "S_BB": S_BB, "S_psi": S_psi, "lag_time": lag_time, "F_BB_peak": F[peak],
              "cond_steps": steps, "cond_time": steps * dt, "largest_B_cluster_fraction_t": cond["largest_B_cluster_fraction"],
              "f_B_t": cond["f_B"], "dense_B_fraction_t": cond["dense_B_fraction"]}

    marks_path = run_dir / "marks.npz"
    if marks_path.is_file():
        with np.load(marks_path) as m:
            mtypes = m["types"]
            mstep = m["steps"].astype(np.int64)
            f_B_marks = m["f_B"]
        mwin = mstep > (1.0 - window_fraction) * production_steps
        C = mark_autocorrelation(mtypes[mwin], max_lag=max(1, int(0.5 * (np.count_nonzero(mwin) - 1))))
        mframe = int(np.min(np.diff(mstep))) if mstep.size > 1 else 1
        mlag = np.arange(C.size) * mframe * dt
        out["mark_memory_time"] = float(relaxation_times(mlag, C[None, :])[0])
        out["mark_persistence"] = persistence_fraction(mtypes[mwin], TYPE_B)
        arrays.update({"mark_lag_time": mlag, "mark_autocorrelation": C, "mark_steps": mstep, "f_B_marks_t": f_B_marks})
    else:
        out["mark_memory_time"] = float("nan")
        out["mark_persistence"] = float("nan")
    return {"summary": out, "arrays": arrays}


def aggregate_campaign(campaign_dir: Path | str, require_complete: bool = True, window_fraction: float = 0.5) -> dict[str, Path]:
    campaign_dir = Path(campaign_dir)
    manifest = json.loads((campaign_dir / "manifest.json").read_text())
    validate_manifest(manifest, DESIGN_KEYS)
    commit = campaign_commit(manifest)
    settings = manifest["simulation_settings"]
    design = manifest["design"]
    dynamic = design["mode"] == "chromatin"
    outputs_required = ("meta.json", "snapshots.csv", "mode_amplitudes.npz", "condensation.csv", "planned_run.json")
    if dynamic:
        outputs_required += ("marks.npz",)
    if settings.get("save_trajectory"):
        outputs_required += ("trajectory.npz",)
    dt = float(settings["dt"])
    production = int(settings["n_steps"])

    analysis_dir = campaign_dir / "analysis"
    (analysis_dir / "runs").mkdir(parents=True, exist_ok=True)
    rows, missing, arrays = [], [], {}
    for spec in manifest["runs"]:
        run_dir = campaign_dir / "runs" / spec["run_id"]
        if not completion_valid(run_dir, spec, outputs_required, commit):
            missing.append(spec["run_id"])
            continue
        cache = analysis_dir / "runs" / f"{spec['run_id']}.npz"
        summary_path = analysis_dir / "runs" / f"{spec['run_id']}.json"
        if cache.is_file() and summary_path.is_file():
            summary = json.loads(summary_path.read_text())
            if summary.get("_window_fraction") != window_fraction:
                summary = None
        else:
            summary = None
        if summary is None:
            result = analyze_run(run_dir, window_fraction, dt, production)
            summary = {**{k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in result["summary"].items()},
                       "_window_fraction": window_fraction}
            np.savez_compressed(cache, **result["arrays"])
            summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        arrays[spec["run_id"]] = dict(np.load(cache))
        rows.append({**{k: spec[k] for k in RUN_KEYS}, "run_id": spec["run_id"],
                     **{k: (float("nan") if v is None else v) for k, v in summary.items() if not k.startswith("_")}})
    if missing and require_complete:
        raise RuntimeError(f"{len(missing)} runs incomplete; pass --allow-incomplete-analysis to aggregate anyway")
    if not rows:
        raise RuntimeError("no completed runs")

    def write_csv(path: Path, fields: list[str], data: list[dict[str, Any]]) -> None:
        with path.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            w.writeheader()
            for r in data:
                w.writerow({k: ("" if isinstance(r.get(k), float) and not np.isfinite(r[k]) else r.get(k, "")) for k in fields})

    out: dict[str, Path] = {}
    per_run_fields = ["run_id", *RUN_KEYS, *OBSERVABLES]
    out["per_run"] = analysis_dir / "per_run.csv"
    write_csv(out["per_run"], per_run_fields, rows)

    cond_keys = ("eps_BB", "kappa", "f_A", "k_off", "k_fb")
    conditions = []
    for key in sorted({tuple(r[k] for k in cond_keys) for r in rows}):
        group = [r for r in rows if tuple(r[k] for k in cond_keys) == key]
        rec = dict(zip(cond_keys, key))
        rec["k_on"] = group[0]["k_on"]
        rec["n_seeds"] = len(group)
        for name in OBSERVABLES:
            m, s = _mean_sem([r[name] for r in group])
            rec[f"{name}_mean"], rec[f"{name}_sem"] = m, s
        rec["f_B_runaway_fraction"] = float(np.mean([(r["f_B_final"] > 0.95) or (r["f_B_final"] < 0.02) for r in group]))
        conditions.append(rec)
    cond_fields = [*cond_keys, "k_on", "n_seeds", "f_B_runaway_fraction"] + [f"{n}_{s}" for n in OBSERVABLES for s in ("mean", "sem")]
    out["per_condition"] = analysis_dir / "per_condition.csv"
    write_csv(out["per_condition"], cond_fields, conditions)
    (analysis_dir / "summary.json").write_text(json.dumps({
        "campaign_id": manifest["campaign_id"], "mode": design["mode"], "n_runs_analyzed": len(rows), "missing_runs": missing,
        "window_fraction": window_fraction, "box_size": manifest["runs"][0]["box_size"], "bead_density": design["bead_density"],
    }, indent=2) + "\n")
    out["summary"] = analysis_dir / "summary.json"
    out.update(make_figures(analysis_dir, conditions, rows, arrays, design, dt))
    return out


def make_figures(analysis_dir: Path, conditions: list[dict], rows: list[dict], arrays: dict, design: dict, dt: float) -> dict[str, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = analysis_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    out: dict[str, Path] = {}
    dynamic = design["mode"] == "chromatin"
    xkey, ykey = ("k_off", "k_fb") if dynamic else ("eps_BB", "kappa")
    xs = sorted({c[xkey] for c in conditions})
    ys = sorted({c[ykey] for c in conditions})
    slices = sorted({tuple(c[k] for k in ("eps_BB", "f_A")) for c in conditions}) if dynamic else sorted({c["f_A"] for c in conditions})

    def heat(metric, label, name, log=False):
        n = len(slices)
        fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 4.2), squeeze=False)
        for ax, sl in zip(axes[0], slices):
            grid = np.full((len(ys), len(xs)), np.nan)
            for c in conditions:
                ok = (c["eps_BB"], c["f_A"]) == sl if dynamic else c["f_A"] == sl
                if ok:
                    grid[ys.index(c[ykey]), xs.index(c[xkey])] = c[metric if metric in c else f"{metric}_mean"]
            im = ax.imshow(grid, origin="lower", aspect="auto", cmap="viridis",
                           norm=matplotlib.colors.LogNorm() if log and np.nanmin(grid) > 0 else None)
            ax.set_xticks(range(len(xs)), [f"{x:g}" for x in xs], rotation=45)
            ax.set_yticks(range(len(ys)), [f"{y:g}" for y in ys])
            ax.set_xlabel(xkey.replace("_", " ")); ax.set_ylabel(ykey.replace("_", " "))
            ax.set_title((f"eps_BB={sl[0]:g}, f_A={sl[1]:g}" if dynamic else f"f_A={sl:g}"), fontsize=9)
            for (i, j), v in np.ndenumerate(grid):
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="w" if (not log and v < np.nanmax(grid) * 0.6) else "k")
            fig.colorbar(im, ax=ax, label=label)
        fig.tight_layout(); fig.savefig(fig_dir / f"{name}.png", dpi=160); fig.savefig(fig_dir / f"{name}.pdf"); plt.close(fig)
        out[f"fig_{name}"] = fig_dir / f"{name}.png"

    heat("largest_B_cluster_fraction", "fraction of B in the largest cluster", "heat_largest_cluster")
    heat("dense_B_fraction", "fraction of B beads with >= n_dense B neighbors", "heat_dense_fraction")
    heat("tau_BB_peak", r"$\tau_{1/e}$ of $F_{BB}(q^*,t)$ [$\tau$]", "heat_blob_lifetime", log=True)
    heat("S_BB_peak", r"$S_{BB}(q^*)$", "heat_S_BB_peak", log=True)
    heat("mean_Rg", r"mean $R_g$", "heat_Rg")
    if dynamic:
        heat("f_B", "mean B fraction (window)", "heat_f_B")
        heat("mark_memory_time", r"mark memory time $\tau_{1/e}$ of $C_m$ [$\tau$]", "heat_mark_memory", log=True)
        heat("f_B_runaway_fraction", "fraction of seeds that ran away (all B or all A)", "heat_runaway")

    # Time series: one panel per slice; curves colored by the x axis, line style by the y axis.
    for sl in slices:
        sel = [r for r in rows if ((r["eps_BB"], r["f_A"]) == sl if dynamic else r["f_A"] == sl)]
        fig, axes = plt.subplots(1, 2 if dynamic else 1, figsize=(11 if dynamic else 6, 4.2), squeeze=False)
        cmap = plt.get_cmap("viridis")
        styles = ["-", "--", ":", "-."]
        seen = set()
        for r in sel:
            if r["seed"] != min(rr["seed"] for rr in sel if rr[xkey] == r[xkey] and rr[ykey] == r[ykey]):
                continue
            a = arrays[r["run_id"]]
            color = cmap((xs.index(r[xkey]) + 0.5) / len(xs)); ls = styles[ys.index(r[ykey]) % len(styles)]
            label = f"{xkey}={r[xkey]:g}, {ykey}={r[ykey]:g}"
            axes[0, 0].plot(a["cond_time"], a["largest_B_cluster_fraction_t"], color=color, ls=ls, lw=1.0, label=label if label not in seen else None)
            if dynamic:
                axes[0, 1].plot(a["mark_steps"] * dt, a["f_B_marks_t"], color=color, ls=ls, lw=1.0)
            seen.add(label)
        axes[0, 0].set_xlabel(r"time ($\tau$)"); axes[0, 0].set_ylabel("largest B cluster fraction"); axes[0, 0].set_xscale("log")
        axes[0, 0].legend(fontsize=6, ncol=2)
        if dynamic:
            axes[0, 1].set_xlabel(r"time ($\tau$)"); axes[0, 1].set_ylabel("B fraction"); axes[0, 1].set_xscale("log")
        title = f"eps_BB={sl[0]:g}_fA={sl[1]:g}" if dynamic else f"fA={sl:g}"
        fig.suptitle(title, fontsize=9); fig.tight_layout()
        fig.savefig(fig_dir / f"timeseries_{title}.png", dpi=160); plt.close(fig)
        out[f"fig_timeseries_{title}"] = fig_dir / f"timeseries_{title}.png"
    return out
