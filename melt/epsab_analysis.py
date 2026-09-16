"""Aggregate an ``eps_AB`` sweep: fluctuation observables versus incompatibility, transition
estimates with seed bootstrap, the RPA comparison, and figures.

The incompatibility axis is ``delta_eps = eps_like - eps_AB`` (zero at the chi = 0 reference).
Per-run observables come from :func:`melt.dynamic_structure.analyze_run`; per-condition
statistics are mean and standard error over seeds.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from .campaign import campaign_commit, completion_valid, validate_manifest
from .dynamic_structure import analyze_run
from .rpa import (
    fit_alpha_from_mixed_side,
    rpa_spinodal_chi,
    rpa_structure_factor,
    segment_length_from_rg,
)

DESIGN_KEYS = ("schema_version", "study", "campaign_id", "design", "simulation_settings", "runs")
RUN_FIELDS = (
    "run_id", "n_chains", "chain_length", "box_size", "eps_AB", "delta_eps", "seed",
    "q_peak", "S_psi_peak", "S_psi_peak_sem_time", "S_psi_kmin", "S_rho_kmin",
    "non_gaussian_ratio_peak", "peak_intensity_relative_variance",
    "cg_variance_1p5", "cg_variance_2", "cg_variance_3",
    "tau_psi_peak", "gamma_psi_peak", "kww_tau_peak", "kww_beta_peak", "plateau_psi_peak",
    "S_psi_peak_relative_drift", "energy_per_bead", "mean_Rg", "T_inst_star", "n_frames_window",
)
CONDITION_OBSERVABLES = (
    "S_psi_peak", "S_psi_kmin", "S_rho_kmin", "non_gaussian_ratio_peak",
    "peak_intensity_relative_variance", "cg_variance_1p5", "cg_variance_2", "cg_variance_3",
    "tau_psi_peak", "gamma_psi_peak", "kww_tau_peak", "kww_beta_peak", "plateau_psi_peak",
    "S_psi_peak_relative_drift", "energy_per_bead", "mean_Rg", "q_peak",
)


def _fmt(value: Any) -> str:
    if isinstance(value, (float, np.floating)):
        x = float(value)
        return "" if not np.isfinite(x) else format(x, ".10g")
    return str(value)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row.get(key, "")) for key in fieldnames})


def _mean_sem(values: np.ndarray) -> tuple[float, float]:
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan")
    sem = float(np.std(x, ddof=1) / np.sqrt(x.size)) if x.size > 1 else float("nan")
    return float(np.mean(x)), sem


def snapshot_window_means(run_dir: Path, first_step: int) -> dict[str, float]:
    energies, rgs, temps = [], [], []
    with (run_dir / "snapshots.csv").open() as handle:
        for row in csv.DictReader(handle):
            if int(row["step"]) < first_step:
                continue
            energies.append(float(row["total_energy"]))
            rgs.append(float(row["mean_Rg"]))
            temps.append(float(row["temperature_inst"]))
    n_beads = None
    with (run_dir / "meta.json").open() as handle:
        n_beads = int(json.load(handle)["n_particles"])
    return {
        "energy_per_bead": float(np.mean(energies)) / n_beads if energies else float("nan"),
        "mean_Rg": float(np.mean(rgs)) if rgs else float("nan"),
        "T_inst_star": float(np.mean(temps)) if temps else float("nan"),
    }


def per_run_row(spec: dict[str, Any], summary: dict[str, Any], window: dict[str, float], eps_like: float) -> dict[str, Any]:
    cg = summary["coarse_grained_variance"]
    return {
        "run_id": spec["run_id"],
        "n_chains": spec["n_chains"],
        "chain_length": spec["chain_length"],
        "box_size": spec["box_size"],
        "eps_AB": spec["eps_AB"],
        "delta_eps": eps_like - float(spec["eps_AB"]),
        "seed": spec["seed"],
        "q_peak": summary["q_peak"],
        "S_psi_peak": summary["S_psi_peak"],
        "S_psi_peak_sem_time": summary["S_psi_peak_sem_time"],
        "S_psi_kmin": summary["S_psi_kmin"],
        "S_rho_kmin": summary["S_rho_kmin"],
        "non_gaussian_ratio_peak": summary["non_gaussian_ratio_peak"],
        "peak_intensity_relative_variance": summary["peak_intensity_relative_variance"],
        "cg_variance_1p5": cg.get("1.5", float("nan")),
        "cg_variance_2": cg.get("2", float("nan")),
        "cg_variance_3": cg.get("3", float("nan")),
        "tau_psi_peak": summary["tau_psi_peak"],
        "gamma_psi_peak": summary["gamma_psi_peak"],
        "kww_tau_peak": summary["kww_psi_peak"]["tau"],
        "kww_beta_peak": summary["kww_psi_peak"]["beta"],
        "plateau_psi_peak": summary["plateau_psi_peak"],
        "S_psi_peak_relative_drift": summary["S_psi_peak_relative_drift"],
        "n_frames_window": summary["n_frames_window"],
        **window,
    }


def condition_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = sorted({(int(r["n_chains"]), float(r["eps_AB"])) for r in rows}, key=lambda k: (k[0], -k[1]))
    out = []
    for n_chains, eps in keys:
        group = [r for r in rows if int(r["n_chains"]) == n_chains and float(r["eps_AB"]) == eps]
        record: dict[str, Any] = {
            "n_chains": n_chains,
            "eps_AB": eps,
            "delta_eps": group[0]["delta_eps"],
            "n_seeds": len(group),
            "box_size": group[0]["box_size"],
        }
        for name in CONDITION_OBSERVABLES:
            mean, sem = _mean_sem(np.asarray([r[name] for r in group], dtype=np.float64))
            record[f"{name}_mean"] = mean
            record[f"{name}_sem"] = sem
        peaks = np.asarray([r["S_psi_peak"] for r in group], dtype=np.float64)
        record["S_psi_peak_seed_relative_variance"] = (
            float(np.var(peaks, ddof=1) / np.mean(peaks) ** 2) if peaks.size > 1 and np.mean(peaks) > 0 else float("nan")
        )
        out.append(record)
    return out


def _argmax_location(x: np.ndarray, y: np.ndarray) -> float:
    good = np.isfinite(y)
    if not np.any(good):
        return float("nan")
    return float(x[good][int(np.argmax(y[good]))])


def _max_log_slope_location(x: np.ndarray, y: np.ndarray) -> float:
    """Midpoint of the interval where ``d ln y / dx`` is largest (x ascending)."""
    good = np.isfinite(y) & (y > 0)
    if np.count_nonzero(good) < 3:
        return float("nan")
    xs, ys = x[good], np.log(y[good])
    slopes = np.diff(ys) / np.diff(xs)
    i = int(np.argmax(slopes))
    return float(0.5 * (xs[i] + xs[i + 1]))


def _crossing_location(x: np.ndarray, y: np.ndarray, level: float) -> float:
    """First x (ascending) where y crosses ``level`` downward, linearly interpolated."""
    good = np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return float("nan")
    xs, ys = x[good], y[good]
    for i in range(xs.size - 1):
        if (ys[i] - level) * (ys[i + 1] - level) <= 0 and ys[i] != ys[i + 1]:
            return float(xs[i] + (level - ys[i]) * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i]))
    return float("nan")


def transition_estimates(rows: list[dict[str, Any]], n_chains: int, n_bootstrap: int = 400, seed: int = 0) -> dict[str, Any]:
    """Locate the transition along ``delta_eps`` with several estimators and a seed bootstrap."""
    rng = np.random.default_rng(seed)
    group = [r for r in rows if int(r["n_chains"]) == n_chains]
    deltas = np.asarray(sorted({float(r["delta_eps"]) for r in group}))
    by_delta = {d: [r for r in group if float(r["delta_eps"]) == d] for d in deltas}

    def curves(sample: bool) -> dict[str, np.ndarray]:
        out = {k: np.full(deltas.size, np.nan) for k in ("S_peak", "rel_var", "R4", "seed_var", "cg2")}
        for i, d in enumerate(deltas):
            rs = by_delta[d]
            if sample:
                rs = [rs[j] for j in rng.integers(0, len(rs), size=len(rs))]
            peaks = np.asarray([r["S_psi_peak"] for r in rs], dtype=np.float64)
            out["S_peak"][i] = np.nanmean(peaks)
            out["rel_var"][i] = np.nanmean([r["peak_intensity_relative_variance"] for r in rs])
            out["R4"][i] = np.nanmean([r["non_gaussian_ratio_peak"] for r in rs])
            out["cg2"][i] = np.nanmean([r["cg_variance_2"] for r in rs])
            if peaks.size > 1 and np.nanmean(peaks) > 0:
                out["seed_var"][i] = np.nanvar(peaks, ddof=1) / np.nanmean(peaks) ** 2
        return out

    def estimators(c: dict[str, np.ndarray]) -> dict[str, float]:
        return {
            "delta_eps_max_peak_intensity_variance": _argmax_location(deltas, c["rel_var"]),
            "delta_eps_max_seed_variance": _argmax_location(deltas, c["seed_var"]),
            "delta_eps_max_log_slope_S_peak": _max_log_slope_location(deltas, c["S_peak"]),
            "delta_eps_max_log_slope_cg_variance": _max_log_slope_location(deltas, c["cg2"]),
            "delta_eps_non_gaussian_ratio_1p5": _crossing_location(deltas, c["R4"], 1.5),
        }

    point = estimators(curves(sample=False))
    boots = {k: [] for k in point}
    for _ in range(int(n_bootstrap)):
        for k, v in estimators(curves(sample=True)).items():
            boots[k].append(v)
    result: dict[str, Any] = {"n_chains": n_chains, "delta_eps_grid": deltas.tolist(), "estimates": {}}
    for k, v in point.items():
        b = np.asarray(boots[k], dtype=np.float64)
        b = b[np.isfinite(b)]
        result["estimates"][k] = {
            "value": v,
            "bootstrap_p16": float(np.percentile(b, 16)) if b.size else float("nan"),
            "bootstrap_p50": float(np.percentile(b, 50)) if b.size else float("nan"),
            "bootstrap_p84": float(np.percentile(b, 84)) if b.size else float("nan"),
            "n_bootstrap_finite": int(b.size),
        }
    return result


def rpa_comparison(
    conditions: list[dict[str, Any]],
    n_chains: int,
    settings: dict[str, Any],
    design: dict[str, Any],
    n_fit_points: int = 4,
) -> dict[str, Any]:
    """Fit the eps_AB -> chi bridge on the mixed side and predict the mean-field spinodal."""
    group = sorted([c for c in conditions if int(c["n_chains"]) == n_chains], key=lambda c: -float(c["eps_AB"]))
    if not group:
        raise ValueError("no conditions for RPA comparison")
    eps_like = float(settings["lj_eps_AA"])
    T_star = float(settings["T_quench"])
    N = int(design["chain_length"])
    kappa = float(design["kappa"])
    pi = float(design["pi"])
    reference = group[0]
    rg = float(reference["mean_Rg_mean"])
    b = segment_length_from_rg(rg, N) if np.isfinite(rg) and rg > 0 else 1.0
    L = float(reference["box_size"])
    q_min = 2.0 * np.pi / L
    q_peak = float(reference["q_peak_mean"]) if np.isfinite(reference["q_peak_mean"]) else q_min
    mixed = [c for c in group if np.isfinite(c["S_psi_peak_mean"]) and c["S_psi_peak_mean"] > 0][:n_fit_points]
    out: dict[str, Any] = {
        "n_chains": n_chains, "eps_like": eps_like, "T_star": T_star, "chain_length": N, "kappa": kappa, "pi": pi,
        "segment_length_b": float(b), "reference_Rg": rg, "q_min": float(q_min), "q_fit": float(q_peak),
        "fit_eps_ABs": [float(c["eps_AB"]) for c in mixed],
    }
    chi_s, q_s = rpa_spinodal_chi(N, kappa, pi, 0.5, b, q=np.asarray([q_min]))
    out["chi_spinodal_at_q_min"] = float(chi_s)
    out["chi_spinodal_q0"] = float(rpa_spinodal_chi(N, kappa, pi, 0.5, b)[0])
    try:
        fit = fit_alpha_from_mixed_side(
            np.asarray([c["eps_AB"] for c in mixed], dtype=np.float64),
            np.asarray([c["S_psi_peak_mean"] for c in mixed], dtype=np.float64),
            N, kappa, pi, 0.5, T_star, b, q_peak, eps_like=eps_like, q_accessible=np.asarray([q_min]),
        )
        out.update({
            "alpha": float(fit["alpha"]),
            "residual_rms_log": float(fit["residual_rms"]),
            "eps_AB_spinodal_rpa": float(fit["eps_AB_spinodal"]),
            "delta_eps_spinodal_rpa": eps_like - float(fit["eps_AB_spinodal"]),
            "n_fit_points": int(fit["n_points"]),
        })
        eps_all = np.asarray([c["eps_AB"] for c in group], dtype=np.float64)
        chi_all = fit["alpha"] * (eps_like - eps_all) / T_star
        out["predicted_S_peak"] = {
            f"{e:g}": (float(v) if np.isfinite(v) else None)
            for e, v in zip(eps_all, rpa_structure_factor(q_peak, chi_all, N, kappa, pi, 0.5, b))
        }
    except (ValueError, RuntimeError) as exc:
        out["fit_error"] = str(exc)
    return out


def _load_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


def aggregate_campaign(campaign_dir: Path | str, require_complete: bool = True, window_fraction: float = 0.5) -> dict[str, Path]:
    campaign_dir = Path(campaign_dir)
    with (campaign_dir / "manifest.json").open() as handle:
        manifest = json.load(handle)
    validate_manifest(manifest, DESIGN_KEYS)
    commit = campaign_commit(manifest)
    settings = manifest["simulation_settings"]
    design = manifest["design"]
    eps_like = float(settings["lj_eps_AA"])
    required = ("meta.json", "snapshots.csv", "structure_factor.npz", "mode_amplitudes.npz", "planned_run.json")
    if settings.get("save_trajectory"):
        required = (*required, "trajectory.npz")

    analysis_dir = campaign_dir / "analysis"
    (analysis_dir / "runs").mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    per_run_arrays: dict[str, dict[str, np.ndarray]] = {}
    for spec in manifest["runs"]:
        run_dir = campaign_dir / "runs" / spec["run_id"]
        if not completion_valid(run_dir, spec, required, commit):
            missing.append(spec["run_id"])
            continue
        cache = analysis_dir / "runs" / spec["run_id"]
        cache.mkdir(parents=True, exist_ok=True)
        summary_path = cache / "fluctuation_summary.json"
        arrays_path = cache / "dynamic_structure.npz"
        if summary_path.is_file() and arrays_path.is_file():
            with summary_path.open() as handle:
                summary = json.load(handle)
            if summary.get("window_fraction") != float(window_fraction):
                summary = None
        else:
            summary = None
        if summary is None:
            result = analyze_run(str(run_dir), window_fraction=window_fraction, write=False)
            summary = _nan_to_none(result["summary"])
            np.savez_compressed(arrays_path, **result["arrays"])
            with summary_path.open("w") as handle:
                json.dump(summary, handle, indent=2, sort_keys=True)
                handle.write("\n")
        summary = _none_to_nan(summary)
        per_run_arrays[spec["run_id"]] = _load_arrays(arrays_path)
        window = snapshot_window_means(run_dir, int(summary["window_first_step"]))
        rows.append(per_run_row(spec, summary, window, eps_like))
    if missing and require_complete:
        raise RuntimeError(f"{len(missing)} runs are incomplete; rerun them or pass --allow-incomplete-analysis")
    if not rows:
        raise RuntimeError("no completed runs to aggregate")

    conditions = condition_table(rows)
    outputs: dict[str, Path] = {}
    outputs["per_run"] = analysis_dir / "per_run.csv"
    _write_csv(outputs["per_run"], list(RUN_FIELDS), rows)
    cond_fields = ["n_chains", "eps_AB", "delta_eps", "n_seeds", "box_size", "S_psi_peak_seed_relative_variance"]
    cond_fields += [f"{n}_{s}" for n in CONDITION_OBSERVABLES for s in ("mean", "sem")]
    outputs["per_condition"] = analysis_dir / "per_condition.csv"
    _write_csv(outputs["per_condition"], cond_fields, conditions)

    spectra_rows, f_rows, tau_rows = [], [], []
    sizes = sorted({int(r["n_chains"]) for r in rows})
    for n_chains in sizes:
        for eps in sorted({float(r["eps_AB"]) for r in rows if int(r["n_chains"]) == n_chains}, reverse=True):
            ids = [r["run_id"] for r in rows if int(r["n_chains"]) == n_chains and float(r["eps_AB"]) == eps]
            arr = [per_run_arrays[i] for i in ids]
            shell_q = arr[0]["shell_q"]
            S_stack = np.stack([a["S_psi_shell"] for a in arr])
            S_mean = S_stack.mean(axis=0)
            S_sem = S_stack.std(axis=0, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else np.full_like(S_mean, np.nan)
            for j in range(shell_q.size):
                spectra_rows.append({"n_chains": n_chains, "eps_AB": eps, "delta_eps": eps_like - eps,
                                     "shell_q": float(shell_q[j]), "S_psi_mean": float(S_mean[j]), "S_psi_sem": float(S_sem[j])})
            peak = int(np.argmax(S_mean))
            lags = arr[0]["lag_time"]
            F_stack = np.stack([a["F_psi_qt"][peak] for a in arr])
            F_mean = F_stack.mean(axis=0)
            F_sem = F_stack.std(axis=0, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else np.full_like(F_mean, np.nan)
            for j in range(lags.size):
                f_rows.append({"n_chains": n_chains, "eps_AB": eps, "delta_eps": eps_like - eps, "q_peak": float(shell_q[peak]),
                               "lag_time": float(lags[j]), "F_psi_mean": float(F_mean[j]), "F_psi_sem": float(F_sem[j])})
            tau_stack = np.stack([a["tau_psi"] for a in arr])
            for j in range(shell_q.size):
                m, s = _mean_sem(tau_stack[:, j])
                tau_rows.append({"n_chains": n_chains, "eps_AB": eps, "delta_eps": eps_like - eps,
                                 "shell_q": float(shell_q[j]), "tau_psi_mean": m, "tau_psi_sem": s})
    outputs["spectra"] = analysis_dir / "spectra_by_condition.csv"
    _write_csv(outputs["spectra"], ["n_chains", "eps_AB", "delta_eps", "shell_q", "S_psi_mean", "S_psi_sem"], spectra_rows)
    outputs["F_peak"] = analysis_dir / "F_peak_by_condition.csv"
    _write_csv(outputs["F_peak"], ["n_chains", "eps_AB", "delta_eps", "q_peak", "lag_time", "F_psi_mean", "F_psi_sem"], f_rows)
    outputs["tau_by_shell"] = analysis_dir / "tau_by_shell_condition.csv"
    _write_csv(outputs["tau_by_shell"], ["n_chains", "eps_AB", "delta_eps", "shell_q", "tau_psi_mean", "tau_psi_sem"], tau_rows)

    summary = {
        "campaign_id": manifest["campaign_id"],
        "design_sha256": manifest["design_sha256"],
        "n_runs_planned": len(manifest["runs"]),
        "n_runs_analyzed": len(rows),
        "missing_runs": missing,
        "window_fraction": float(window_fraction),
        "incompatibility_axis": f"delta_eps = {eps_like:g} - eps_AB",
        "transition": {str(n): transition_estimates(rows, n) for n in sizes},
        "rpa": {str(n): rpa_comparison(conditions, n, settings, design) for n in sizes},
    }
    outputs["transition_summary"] = analysis_dir / "transition_summary.json"
    with outputs["transition_summary"].open("w") as handle:
        json.dump(_nan_to_none(summary), handle, indent=2, sort_keys=True)
        handle.write("\n")
    outputs.update(make_figures(analysis_dir, conditions, spectra_rows, f_rows, summary))
    return outputs


def _nan_to_none(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _nan_to_none(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_nan_to_none(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def _none_to_nan(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _none_to_nan(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_none_to_nan(v) for v in value]
    return float("nan") if value is None else value


def make_figures(analysis_dir: Path, conditions, spectra_rows, f_rows, summary) -> dict[str, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = analysis_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    sizes = sorted({int(c["n_chains"]) for c in conditions})

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(fig_dir / f"{name}.png", dpi=160)
        fig.savefig(fig_dir / f"{name}.pdf")
        plt.close(fig)
        outputs[f"fig_{name}"] = fig_dir / f"{name}.png"

    def series(n_chains, key):
        g = sorted([c for c in conditions if int(c["n_chains"]) == n_chains], key=lambda c: c["delta_eps"])
        x = np.asarray([c["delta_eps"] for c in g], dtype=np.float64)
        y = np.asarray([c[f"{key}_mean"] for c in g], dtype=np.float64)
        e = np.asarray([c[f"{key}_sem"] for c in g], dtype=np.float64)
        return x, y, e, g

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for n in sizes:
        x, y, e, g = series(n, "S_psi_peak")
        ax.errorbar(x, y, yerr=np.nan_to_num(e), marker="o", ms=4, lw=1.2, capsize=2, label=f"M={n} simulation")
        rpa = summary["rpa"].get(str(n), {})
        pred = rpa.get("predicted_S_peak")
        if pred:
            eps_like = rpa["eps_like"]
            xs = np.asarray([eps_like - float(k) for k in pred])
            ys = np.asarray([np.nan if v is None else v for v in pred.values()], dtype=np.float64)
            order = np.argsort(xs)
            ax.plot(xs[order], ys[order], "--", lw=1.0, label=f"RPA fit (alpha={rpa.get('alpha', float('nan')):.2f})")
            if rpa.get("delta_eps_spinodal_rpa") is not None:
                ax.axvline(rpa["delta_eps_spinodal_rpa"], color="0.5", ls=":", lw=1.0, label="RPA spinodal")
        est = summary["transition"][str(n)]["estimates"]
        v = est["delta_eps_max_peak_intensity_variance"]["value"]
        if v is not None and np.isfinite(v):
            ax.axvline(v, color="C3", ls="-.", lw=1.0, label="max intensity variance")
    ax.set_yscale("log")
    ax.set_xlabel(r"incompatibility $\Delta\varepsilon=\varepsilon_{AA}-\varepsilon_{AB}$")
    ax.set_ylabel(r"$S_{\psi\psi}(q^*)$ (time-averaged, window)")
    ax.legend(fontsize=8)
    save(fig, "S_peak_vs_delta_eps")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
    for n in sizes:
        for ax, key, label in zip(
            axes,
            ("peak_intensity_relative_variance", "cg_variance_2", "non_gaussian_ratio_peak"),
            (r"Var$_t[S(q^*,t)]/\langle S\rangle^2$", r"coarse-grained variance $\ell=2\sigma$", r"$\langle|\rho_\psi|^4\rangle/\langle|\rho_\psi|^2\rangle^2$"),
        ):
            x, y, e, _ = series(n, key)
            ax.errorbar(x, y, yerr=np.nan_to_num(e), marker="o", ms=4, lw=1.2, capsize=2, label=f"M={n}")
            ax.set_xlabel(r"$\Delta\varepsilon$")
            ax.set_ylabel(label)
        g = sorted([c for c in conditions if int(c["n_chains"]) == n], key=lambda c: c["delta_eps"])
        axes[0].plot([c["delta_eps"] for c in g], [c["S_psi_peak_seed_relative_variance"] for c in g], "s--", ms=3, lw=0.8, label=f"M={n} seed variance")
    axes[2].axhline(2.0, color="0.6", lw=0.8, ls=":")
    axes[2].axhline(1.0, color="0.6", lw=0.8, ls=":")
    axes[1].set_yscale("log")
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=7)
    save(fig, "fluctuation_measures_vs_delta_eps")

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    for n in sizes:
        x, y, e, _ = series(n, "tau_psi_peak")
        axes[0].errorbar(x, y, yerr=np.nan_to_num(e), marker="o", ms=4, lw=1.2, capsize=2, label=f"M={n}")
        x, y, e, _ = series(n, "plateau_psi_peak")
        axes[1].errorbar(x, y, yerr=np.nan_to_num(e), marker="o", ms=4, lw=1.2, capsize=2, label=f"M={n}")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$\Delta\varepsilon$")
    axes[0].set_ylabel(r"$\tau_{1/e}(q^*)$ of $F_{\psi\psi}(q^*,t)$")
    axes[1].set_xlabel(r"$\Delta\varepsilon$")
    axes[1].set_ylabel(r"late-lag plateau of $F_{\psi\psi}(q^*,t)$")
    axes[0].legend(fontsize=8)
    save(fig, "relaxation_vs_delta_eps")

    for n in sizes:
        deltas = sorted({r["delta_eps"] for r in f_rows if int(r["n_chains"]) == n})
        pick = deltas if len(deltas) <= 10 else [deltas[i] for i in np.linspace(0, len(deltas) - 1, 10).astype(int)]
        cmap = plt.get_cmap("viridis")
        fig, ax = plt.subplots(figsize=(6.0, 4.2))
        for d in pick:
            rs = [r for r in f_rows if int(r["n_chains"]) == n and r["delta_eps"] == d]
            t = np.asarray([r["lag_time"] for r in rs]); F = np.asarray([r["F_psi_mean"] for r in rs])
            color = cmap((d - deltas[0]) / max(deltas[-1] - deltas[0], 1e-12))
            ax.plot(t[1:], F[1:], color=color, lw=1.2, label=rf"$\Delta\varepsilon$={d:g} (q*={rs[0]['q_peak']:.2f})")
        ax.set_xscale("log")
        ax.set_xlabel(r"lag $t$ ($\tau$)")
        ax.set_ylabel(r"$F_{\psi\psi}(q^*,t)=S(q^*,t)/S(q^*,0)$")
        ax.set_ylim(-0.1, 1.05)
        ax.axhline(0.0, color="0.7", lw=0.8)
        ax.legend(fontsize=6, ncol=2)
        save(fig, f"F_peak_vs_lag_M{n}")

        fig, ax = plt.subplots(figsize=(6.0, 4.2))
        for d in pick:
            rs = [r for r in spectra_rows if int(r["n_chains"]) == n and r["delta_eps"] == d]
            q = np.asarray([r["shell_q"] for r in rs]); S = np.asarray([r["S_psi_mean"] for r in rs])
            color = cmap((d - deltas[0]) / max(deltas[-1] - deltas[0], 1e-12))
            ax.plot(q, S, marker="o", ms=3, lw=1.0, color=color, label=rf"$\Delta\varepsilon$={d:g}")
        ax.set_yscale("log")
        ax.set_xlabel(r"$q$ ($\sigma^{-1}$)")
        ax.set_ylabel(r"$S_{\psi\psi}(q)$")
        ax.legend(fontsize=6, ncol=2)
        save(fig, f"spectra_M{n}")
    return outputs
