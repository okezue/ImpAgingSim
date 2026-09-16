"""Static fluctuation observables and the dynamic structure factor from recorded box modes.

Input is a ``mode_amplitudes.npz`` written by :class:`melt.modes.ModeAmplitudeRecorder`:
``rho_A(q,t)`` and ``rho_B(q,t)`` at every reciprocal box mode with ``|q|<=q_max``.

Composition channel ``psi = rho_A - rho_B`` (the label-symmetric order parameter) and total
density channel ``rho = rho_A + rho_B``.  All spectra use the per-bead normalization
``S(q)=|rho(q)|^2/N_total`` of :mod:`melt.direct_structure`.  The two-time object is

    S(q,tau) = < rho(q,t0) rho*(q,t0+tau) > / N_total,   F(q,tau) = S(q,tau)/S(q,0),

averaged over time origins ``t0`` inside a stationary window and over the modes of each
cubic shell ``|h|^2 = const``.  Time means are *not* subtracted before correlating, so a
frozen microphase pattern shows up as a nonzero plateau ``F(q,tau->inf)`` instead of being
silently removed.
"""
from __future__ import annotations

import json
import os
from typing import Any

import numpy as np

from .dynamics import kww_fit, tau_alpha_from_overlap
from .modes import load_mode_amplitudes

DEFAULT_COARSE_GRAIN_LENGTHS = (1.5, 2.0, 3.0)


def composition_amplitudes(data: dict[str, np.ndarray]) -> np.ndarray:
    return np.asarray(data["rho_A"], dtype=np.complex128) - np.asarray(data["rho_B"], dtype=np.complex128)


def density_amplitudes(data: dict[str, np.ndarray]) -> np.ndarray:
    return np.asarray(data["rho_A"], dtype=np.complex128) + np.asarray(data["rho_B"], dtype=np.complex128)


def stationary_window_mask(steps: np.ndarray, production_steps: int, fraction: float) -> np.ndarray:
    """Frames with ``step > (1-fraction)*production_steps``: the trailing part of the run."""
    s = np.asarray(steps, dtype=np.int64)
    if not 0.0 < float(fraction) <= 1.0:
        raise ValueError("window fraction must lie in (0,1]")
    start = (1.0 - float(fraction)) * float(production_steps)
    mask = s > start + 1e-9
    if np.count_nonzero(mask) < 2:
        raise ValueError("stationary window must contain at least two frames")
    return mask


def equal_time_spectrum(rho: np.ndarray, n_total: int) -> np.ndarray:
    return (np.abs(np.asarray(rho, dtype=np.complex128)) ** 2) / float(n_total)


def shell_average_series(values: np.ndarray, shell_index: np.ndarray, n_shells: int) -> np.ndarray:
    """Average a (T,Q) mode-resolved series within shells, returning (T,S)."""
    x = np.asarray(values, dtype=np.float64)
    idx = np.asarray(shell_index, dtype=np.int64)
    if x.ndim != 2 or idx.shape != (x.shape[1],):
        raise ValueError("values must be (T,Q) with one shell index per mode")
    counts = np.bincount(idx, minlength=int(n_shells)).astype(np.float64)
    if counts.size != int(n_shells) or np.any(counts == 0):
        raise ValueError("shell_index must populate every shell")
    sums = np.zeros((x.shape[0], int(n_shells)), dtype=np.float64)
    np.add.at(sums.T, idx, x.T)
    return sums / counts


def block_sem(series: np.ndarray, n_blocks: int = 5) -> np.ndarray:
    """Standard error of the time mean from block averages along axis 0."""
    x = np.asarray(series, dtype=np.float64)
    T = x.shape[0]
    blocks = int(min(n_blocks, T))
    if blocks < 2:
        return np.full(x.shape[1:], np.nan)
    edges = np.linspace(0, T, blocks + 1).astype(int)
    means = np.stack([x[a:b].mean(axis=0) for a, b in zip(edges[:-1], edges[1:])], axis=0)
    return np.std(means, axis=0, ddof=1) / np.sqrt(blocks)


def coarse_grained_variance(
    S_mean_modes: np.ndarray, q: np.ndarray, n_total: int, ell: float
) -> float:
    """Relative variance of the Gaussian-filtered composition field at length ``ell``.

    With ``psi_ell(q)=rho_psi(q)exp(-q^2 ell^2/2)`` the spatial variance divided by the
    squared mean bead density is ``(1/N) sum_{q!=0} S(q) exp(-q^2 ell^2)``.  The sum runs
    over the recorded modes only, so ``ell`` must be large enough that the ``q_max`` cutoff
    is negligible (``exp(-q_max^2 ell^2) << 1``).
    """
    S = np.asarray(S_mean_modes, dtype=np.float64)
    qq = np.asarray(q, dtype=np.float64)
    if S.shape != qq.shape:
        raise ValueError("S_mean_modes and q must align")
    return float(np.sum(S * np.exp(-(qq**2) * float(ell) ** 2)) / float(n_total))


def non_gaussian_ratio(rho_modes: np.ndarray) -> float:
    """``<|rho|^4>/<|rho|^2>^2`` pooled over time and modes: 2 for a complex Gaussian, 1 if frozen."""
    a2 = np.abs(np.asarray(rho_modes, dtype=np.complex128)) ** 2
    m2 = float(np.mean(a2))
    if m2 <= 0.0:
        return float("nan")
    return float(np.mean(a2 * a2) / (m2 * m2))


def time_correlation_by_shell(
    rho: np.ndarray,
    shell_index: np.ndarray,
    n_shells: int,
    n_total: int,
    max_lag: int | None = None,
) -> np.ndarray:
    """``S(q_s,tau)`` for every shell and integer frame lag ``tau`` in ``[0,max_lag]``.

    Uses the Wiener-Khinchin identity: the lagged product summed over time origins is the
    inverse FFT of ``|FFT(rho)|^2`` on a zero-padded axis.  The result is divided by the
    number of time-origin pairs per lag, then by ``N_total``, and averaged within shells.
    """
    x = np.asarray(rho, dtype=np.complex128)
    if x.ndim != 2:
        raise ValueError("rho must be (T,Q)")
    T, Q = x.shape
    lag_max = T - 1 if max_lag is None else int(max_lag)
    if lag_max < 0 or lag_max > T - 1:
        raise ValueError(f"max_lag must lie in [0,{T-1}]")
    n_fft = 1 << int(np.ceil(np.log2(2 * T)))
    spectrum = np.fft.fft(x, n=n_fft, axis=0)
    corr = np.fft.ifft(np.abs(spectrum) ** 2, axis=0)[: lag_max + 1].real
    pairs = (T - np.arange(lag_max + 1)).astype(np.float64)
    corr = corr / pairs[:, None] / float(n_total)
    return shell_average_series(corr, shell_index, n_shells).T


def relaxation_times(lag_time: np.ndarray, F: np.ndarray, threshold: float = 1.0 / np.e) -> np.ndarray:
    """First-crossing relaxation time of each row of ``F`` (NaN if never below threshold)."""
    return np.asarray([tau_alpha_from_overlap(lag_time, row, threshold=threshold) for row in F])


def initial_decay_rates(lag_time: np.ndarray, F: np.ndarray, n_points: int = 4) -> np.ndarray:
    """Short-time rate ``Gamma`` from a linear fit of ``-ln F`` through the origin."""
    t = np.asarray(lag_time, dtype=np.float64)[1 : n_points + 1]
    out = np.full(F.shape[0], np.nan)
    for i, row in enumerate(np.asarray(F, dtype=np.float64)):
        y = row[1 : n_points + 1]
        good = np.isfinite(y) & (y > 0.0)
        if np.count_nonzero(good) < 2:
            continue
        tt = t[good]
        out[i] = float(np.dot(tt, -np.log(y[good])) / np.dot(tt, tt))
    return out


def kww_by_shell(lag_time: np.ndarray, F: np.ndarray) -> np.ndarray:
    """Stretched-exponential ``A exp(-(t/tau)^beta)`` fits per shell, shape (S,3)."""
    return np.asarray([kww_fit(lag_time, row) for row in F], dtype=np.float64)


def analyze_mode_amplitudes(
    data: dict[str, np.ndarray],
    window_fraction: float = 0.5,
    max_lag_fraction: float = 0.5,
    coarse_grain_lengths: tuple[float, ...] = DEFAULT_COARSE_GRAIN_LENGTHS,
    n_blocks: int = 5,
) -> dict[str, Any]:
    steps = np.asarray(data["steps"], dtype=np.int64)
    production_steps = int(data["production_steps"]) if "production_steps" in data else int(steps[-1])
    dt = float(data["dt_ps"]) if "dt_ps" in data else 1.0
    n_total = int(data["n_total"])
    shell_index = np.asarray(data["shell_index"], dtype=np.int64)
    shell_q = np.asarray(data["shell_q"], dtype=np.float64)
    q = np.asarray(data["q"], dtype=np.float64)
    n_shells = shell_q.size

    window = stationary_window_mask(steps, production_steps, window_fraction)
    psi = composition_amplitudes(data)[window]
    rho = density_amplitudes(data)[window]
    steps_w = steps[window]
    T = steps_w.size
    frame_interval = int(np.min(np.diff(steps_w))) if T > 1 else 0
    if T > 1 and np.any(np.diff(steps_w) != frame_interval):
        raise ValueError("recorded mode frames are not uniformly spaced inside the window")

    S_psi_modes = equal_time_spectrum(psi, n_total)
    S_rho_modes = equal_time_spectrum(rho, n_total)
    S_psi_shell_t = shell_average_series(S_psi_modes, shell_index, n_shells)
    S_rho_shell_t = shell_average_series(S_rho_modes, shell_index, n_shells)
    S_psi_shell = S_psi_shell_t.mean(axis=0)
    S_rho_shell = S_rho_shell_t.mean(axis=0)
    S_psi_shell_sem = block_sem(S_psi_shell_t, n_blocks)
    S_psi_mode_mean = S_psi_modes.mean(axis=0)

    peak = int(np.argmax(S_psi_shell))
    peak_modes = shell_index == peak
    half = T // 2
    first_half = float(S_psi_shell_t[:half, peak].mean()) if half > 0 else float("nan")
    second_half = float(S_psi_shell_t[half:, peak].mean()) if half > 0 else float("nan")

    max_lag = int(np.floor(max_lag_fraction * (T - 1)))
    S_psi_qt = time_correlation_by_shell(psi, shell_index, n_shells, n_total, max_lag=max_lag)
    S_rho_qt = time_correlation_by_shell(rho, shell_index, n_shells, n_total, max_lag=max_lag)
    lag_steps = np.arange(max_lag + 1, dtype=np.int64) * frame_interval
    lag_time = lag_steps.astype(np.float64) * dt
    with np.errstate(invalid="ignore", divide="ignore"):
        F_psi = S_psi_qt / S_psi_qt[:, :1]
        F_rho = S_rho_qt / S_rho_qt[:, :1]

    tau_psi = relaxation_times(lag_time, F_psi)
    gamma_psi = initial_decay_rates(lag_time, F_psi)
    kww_psi = kww_by_shell(lag_time, F_psi)
    plateau_psi = F_psi[:, max(1, max_lag // 2) :].mean(axis=1) if max_lag >= 1 else np.full(n_shells, np.nan)

    summary: dict[str, Any] = {
        "n_total": n_total,
        "production_steps": production_steps,
        "window_fraction": float(window_fraction),
        "window_first_step": int(steps_w[0]),
        "window_last_step": int(steps_w[-1]),
        "n_frames_window": int(T),
        "frame_interval_steps": frame_interval,
        "dt": dt,
        "q_min": float(shell_q[0]),
        "peak_shell_index": peak,
        "q_peak": float(shell_q[peak]),
        "S_psi_peak": float(S_psi_shell[peak]),
        "S_psi_peak_sem_time": float(S_psi_shell_sem[peak]),
        "S_psi_kmin": float(S_psi_shell[0]),
        "S_rho_kmin": float(S_rho_shell[0]),
        "S_psi_peak_first_half": first_half,
        "S_psi_peak_second_half": second_half,
        "S_psi_peak_relative_drift": (
            float((second_half - first_half) / (0.5 * (first_half + second_half)))
            if np.isfinite(first_half) and np.isfinite(second_half) and (first_half + second_half) > 0
            else float("nan")
        ),
        "peak_intensity_relative_variance": float(
            np.var(S_psi_shell_t[:, peak], ddof=1) / S_psi_shell[peak] ** 2
        )
        if T > 1 and S_psi_shell[peak] > 0
        else float("nan"),
        "non_gaussian_ratio_peak": non_gaussian_ratio(psi[:, peak_modes]),
        "non_gaussian_ratio_kmin": non_gaussian_ratio(psi[:, shell_index == 0]),
        "coarse_grained_variance": {
            f"{float(ell):g}": coarse_grained_variance(S_psi_mode_mean, q, n_total, ell)
            for ell in coarse_grain_lengths
        },
        "tau_psi_peak": float(tau_psi[peak]),
        "tau_psi_kmin": float(tau_psi[0]),
        "gamma_psi_peak": float(gamma_psi[peak]),
        "kww_psi_peak": {"A": float(kww_psi[peak, 0]), "tau": float(kww_psi[peak, 1]), "beta": float(kww_psi[peak, 2])},
        "plateau_psi_peak": float(plateau_psi[peak]),
        "max_lag_steps": int(lag_steps[-1]),
    }
    arrays = {
        "shell_q": shell_q,
        "shell_n2": np.asarray(data["shell_n2"]) if "shell_n2" in data else np.arange(n_shells),
        "S_psi_shell": S_psi_shell,
        "S_psi_shell_sem_time": S_psi_shell_sem,
        "S_rho_shell": S_rho_shell,
        "S_psi_shell_t": S_psi_shell_t,
        "window_steps": steps_w,
        "lag_steps": lag_steps,
        "lag_time": lag_time,
        "S_psi_qt": S_psi_qt,
        "F_psi_qt": F_psi,
        "S_rho_qt": S_rho_qt,
        "F_rho_qt": F_rho,
        "tau_psi": tau_psi,
        "gamma_psi": gamma_psi,
        "kww_psi": kww_psi,
        "plateau_psi": plateau_psi,
    }
    return {"summary": summary, "arrays": arrays}


def analyze_run(
    run_dir: str,
    window_fraction: float = 0.5,
    max_lag_fraction: float = 0.5,
    coarse_grain_lengths: tuple[float, ...] = DEFAULT_COARSE_GRAIN_LENGTHS,
    write: bool = True,
) -> dict[str, Any]:
    """Analyze ``<run_dir>/mode_amplitudes.npz``; optionally write the derived files next to it."""
    path = os.path.join(run_dir, "mode_amplitudes.npz")
    data = load_mode_amplitudes(path)
    result = analyze_mode_amplitudes(
        data,
        window_fraction=window_fraction,
        max_lag_fraction=max_lag_fraction,
        coarse_grain_lengths=coarse_grain_lengths,
    )
    for key in ("kappa", "pi", "f_A", "lj_eps_AB", "T_quench_star", "seed"):
        if key in data:
            result["summary"][key] = float(data[key]) if key != "seed" else int(data[key])
    if write:
        np.savez_compressed(os.path.join(run_dir, "dynamic_structure.npz"), **result["arrays"])
        with open(os.path.join(run_dir, "fluctuation_summary.json"), "w") as handle:
            json.dump(_jsonable(result["summary"]), handle, indent=2, sort_keys=True)
            handle.write("\n")
    return result


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.floating, float)):
        x = float(value)
        return x if np.isfinite(x) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Dynamic structure factor from recorded box modes")
    parser.add_argument("run_dirs", nargs="+")
    parser.add_argument("--window-fraction", type=float, default=0.5)
    parser.add_argument("--max-lag-fraction", type=float, default=0.5)
    args = parser.parse_args(argv)
    for run_dir in args.run_dirs:
        result = analyze_run(run_dir, args.window_fraction, args.max_lag_fraction)
        s = result["summary"]
        print(
            f"{run_dir}: q*={s['q_peak']:.3f} S_psi(q*)={s['S_psi_peak']:.3f} "
            f"R4={s['non_gaussian_ratio_peak']:.3f} tau(q*)={s['tau_psi_peak']:.3g} "
            f"plateau={s['plateau_psi_peak']:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
