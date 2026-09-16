from __future__ import annotations

import json
import os

import numpy as np
import pytest

from melt.direct_structure import direct_partial_structure_factors, reciprocal_modes
from melt.dynamic_structure import (
    analyze_mode_amplitudes,
    analyze_run,
    coarse_grained_variance,
    initial_decay_rates,
    non_gaussian_ratio,
    relaxation_times,
    shell_anisotropy,
    shell_average_series,
    stationary_window_mask,
    time_correlation_by_shell,
)
from melt.modes import (
    ModeAmplitudeRecorder,
    SeparableModeEvaluator,
    load_mode_amplitudes,
    partial_structure_factors_from_amplitudes,
)

try:
    import openmm  # noqa: F401

    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
needs_openmm = pytest.mark.skipif(not HAS_OPENMM, reason="openmm not installed")


class TestSeparableModes:
    def test_matches_direct_estimator(self):
        rng = np.random.default_rng(0)
        L = 9.0
        N = 300
        pos = rng.uniform(-3.0, L + 3.0, size=(N, 3))
        types = (rng.random(N) < 0.4).astype(np.int8)
        modes = reciprocal_modes(L, 2.0)
        ref = direct_partial_structure_factors(pos, types, modes.q_vectors)
        rho_A, rho_B = SeparableModeEvaluator(L, modes).evaluate(pos, types)
        got = partial_structure_factors_from_amplitudes(rho_A, rho_B, N)
        for key in ("S_AA", "S_BB", "S_AB", "S_CC"):
            np.testing.assert_allclose(got[key], ref[key], rtol=1e-10, atol=1e-10)

    def test_zero_mode_of_single_species_counts_beads(self):
        rng = np.random.default_rng(1)
        L = 7.0
        pos = rng.uniform(0.0, L, size=(50, 3))
        types = np.ones(50, dtype=np.int8)
        modes = reciprocal_modes(L, 1.0)
        ev = SeparableModeEvaluator(L, modes)
        rho_A, rho_B = ev.evaluate(pos, types)
        assert np.allclose(rho_B, 0.0)
        # A pure plane wave of beads at x=0 gives rho_A(h_x,0,0)=N exactly.
        pos0 = np.zeros((50, 3))
        rho_A0, _ = ev.evaluate(pos0, types)
        assert np.allclose(rho_A0, 50.0)

    def test_rejects_bad_types(self):
        modes = reciprocal_modes(5.0, 1.5)
        ev = SeparableModeEvaluator(5.0, modes)
        with pytest.raises(ValueError):
            ev.evaluate(np.zeros((3, 3)), np.array([0, 1, 2]))


class TestRecorder:
    def test_records_saves_and_loads(self, tmp_path):
        rng = np.random.default_rng(2)
        L = 6.0
        N = 40
        types = (rng.random(N) < 0.5).astype(np.int8)
        rec = ModeAmplitudeRecorder(L, 1.6)
        for step in (100, 200, 300):
            rec.observe(step, rng.uniform(0, L, size=(N, 3)), types)
        with pytest.raises(ValueError):
            rec.observe(300, rng.uniform(0, L, size=(N, 3)), types)
        path = rec.save(str(tmp_path / "mode_amplitudes.npz"), production_steps=300, dt_ps=0.005)
        data = load_mode_amplitudes(path)
        assert data["steps"].tolist() == [100, 200, 300]
        assert data["rho_A"].shape == (3, data["q"].size)
        assert data["rho_A"].dtype == np.complex64
        assert int(data["n_total"]) == N
        assert int(data["production_steps"]) == 300

    def test_complex64_storage_error_is_small(self):
        rng = np.random.default_rng(3)
        L = 22.0
        N = 5760
        pos = rng.uniform(0, L, size=(N, 3))
        types = (rng.random(N) < 0.5).astype(np.int8)
        rec = ModeAmplitudeRecorder(L, 1.5)
        rec.observe(1, pos, types)
        payload = rec.payload()
        ref = direct_partial_structure_factors(pos, types, rec.modes.q_vectors)
        got = partial_structure_factors_from_amplitudes(payload["rho_A"][0], payload["rho_B"][0], N)
        assert np.max(np.abs(got["S_CC"] - ref["S_CC"])) < 1e-4 * np.max(ref["S_CC"])


def _ou_modes(rng, T, Q, tau_frames, amplitude2):
    """Stationary complex Ornstein-Uhlenbeck modes with <|rho|^2>=amplitude2 and memory tau."""
    a = np.exp(-1.0 / tau_frames)
    noise_scale = np.sqrt(amplitude2 * (1.0 - a * a) / 2.0)
    x = np.empty((T, Q), dtype=np.complex128)
    x[0] = np.sqrt(amplitude2 / 2.0) * (rng.normal(size=Q) + 1j * rng.normal(size=Q))
    for t in range(1, T):
        x[t] = a * x[t - 1] + noise_scale * (rng.normal(size=Q) + 1j * rng.normal(size=Q))
    return x


class TestTimeCorrelation:
    def test_ou_process_recovers_exponential_decay(self):
        rng = np.random.default_rng(4)
        T, Q, tau = 4000, 60, 12.0
        x = _ou_modes(rng, T, Q, tau, amplitude2=5.0)
        shell_index = np.zeros(Q, dtype=np.int64)
        S = time_correlation_by_shell(x, shell_index, 1, n_total=1, max_lag=60)
        F = S / S[:, :1]
        lags = np.arange(61, dtype=float)
        expected = np.exp(-lags / tau)
        assert abs(S[0, 0] - 5.0) < 0.3
        np.testing.assert_allclose(F[0, :30], expected[:30], atol=0.05)
        tau_est = relaxation_times(lags, F)[0]
        assert abs(tau_est - tau) < 1.0
        gamma = initial_decay_rates(lags, F, n_points=3)[0]
        assert abs(gamma - 1.0 / tau) < 0.02

    def test_matches_direct_lagged_average(self):
        rng = np.random.default_rng(5)
        T, Q = 50, 7
        x = rng.normal(size=(T, Q)) + 1j * rng.normal(size=(T, Q))
        shell_index = np.array([0, 0, 1, 1, 1, 2, 2])
        S = time_correlation_by_shell(x, shell_index, 3, n_total=2, max_lag=10)
        for lag in (0, 3, 10):
            direct = np.real(x[: T - lag] * np.conj(x[lag:])).mean(axis=0) / 2.0
            expected = np.array([direct[shell_index == s].mean() for s in range(3)])
            np.testing.assert_allclose(S[:, lag], expected, rtol=1e-10, atol=1e-10)

    def test_frozen_mode_keeps_plateau(self):
        x = np.full((100, 3), 2.0 + 1.0j)
        S = time_correlation_by_shell(x, np.zeros(3, dtype=int), 1, n_total=1, max_lag=20)
        np.testing.assert_allclose(S[0], 5.0)


class TestStaticObservables:
    def test_non_gaussian_ratio_limits(self):
        rng = np.random.default_rng(6)
        gaussian = rng.normal(size=(20000, 4)) + 1j * rng.normal(size=(20000, 4))
        assert abs(non_gaussian_ratio(gaussian) - 2.0) < 0.05
        frozen = np.exp(1j * rng.uniform(0, 2 * np.pi, size=(200, 4)))
        assert abs(non_gaussian_ratio(frozen) - 1.0) < 1e-12

    def test_lamella_is_anisotropic_not_non_gaussian(self):
        # Six-mode shell, all intensity frozen in one +/- pair: ratio stays 1, anisotropy is 3.
        modes = np.zeros((50, 6), dtype=np.complex128)
        modes[:, 0] = 3.0 * np.exp(0.2j)
        modes[:, 1] = 3.0 * np.exp(-0.2j)
        assert abs(non_gaussian_ratio(modes) - 1.0) < 1e-12
        assert abs(shell_anisotropy(modes) - 3.0) < 1e-12
        rng = np.random.default_rng(12)
        iso = rng.normal(size=(20000, 6)) + 1j * rng.normal(size=(20000, 6))
        assert abs(shell_anisotropy(iso) - 1.0) < 0.05

    def test_coarse_grained_variance_reduces_to_parseval_sum(self):
        q = np.array([0.5, 1.0, 1.5])
        S = np.array([3.0, 2.0, 1.0])
        assert abs(coarse_grained_variance(S, q, 10, 0.0) - 0.6) < 1e-12
        assert coarse_grained_variance(S, q, 10, 2.0) < coarse_grained_variance(S, q, 10, 1.0)

    def test_shell_average_series_and_window(self):
        vals = np.array([[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]])
        out = shell_average_series(vals, np.array([0, 0, 1]), 2)
        np.testing.assert_allclose(out, [[2.0, 5.0], [3.0, 6.0]])
        steps = np.arange(100, 1100, 100)
        mask = stationary_window_mask(steps, 1000, 0.5)
        assert steps[mask].tolist() == [600, 700, 800, 900, 1000]


class TestAnalyzeRun:
    def _synthetic_run(self, tmp_path, rng, tau_frames=8.0, frozen=False):
        L = 8.0
        modes = reciprocal_modes(L, 1.6)
        Q = modes.q.size
        T = 400
        # Shell 1 (second-lowest q) carries the peak; everything else is weaker.
        amp = np.where(modes.shell_index == 1, 20.0, 2.0)
        psi = _ou_modes(rng, T, Q, tau_frames, 1.0) * np.sqrt(amp)
        if frozen:
            psi[:, modes.shell_index == 1] = np.sqrt(20.0) * np.exp(1j * 0.3)
        rho_A = 0.5 * psi
        rho_B = -0.5 * psi
        payload = {
            "steps": np.arange(1, T + 1) * 100,
            "rho_A": rho_A.astype(np.complex64),
            "rho_B": rho_B.astype(np.complex64),
            "h": modes.h,
            "q_vectors": modes.q_vectors,
            "q": modes.q,
            "n2": modes.n2,
            "shell_n2": modes.shell_n2,
            "shell_q": modes.shell_q,
            "shell_index": modes.shell_index,
            "shell_degeneracy": modes.shell_degeneracy,
            "n_A": 50,
            "n_B": 50,
            "n_total": 100,
            "box_size_nm": L,
            "production_steps": T * 100,
            "dt_ps": 0.005,
            "lj_eps_AB": 0.5,
            "kappa": 0.5,
            "pi": 0.99,
            "f_A": 0.5,
            "T_quench_star": 0.7,
            "seed": 3,
        }
        run_dir = tmp_path / "run"
        run_dir.mkdir()
        np.savez_compressed(run_dir / "mode_amplitudes.npz", **payload)
        return str(run_dir), modes

    def test_peak_shell_and_relaxation_time(self, tmp_path):
        rng = np.random.default_rng(7)
        run_dir, modes = self._synthetic_run(tmp_path, rng, tau_frames=8.0)
        result = analyze_run(run_dir)
        s = result["summary"]
        assert s["peak_shell_index"] == 1
        assert abs(s["q_peak"] - modes.shell_q[1]) < 1e-12
        assert abs(s["S_psi_peak"] - 20.0 / 100.0) < 0.05
        # lag_time = frames * 100 steps * 0.005 -> tau = 8 frames = 4.0 time units
        assert abs(s["tau_psi_peak"] - 4.0) < 1.0
        assert 1.7 < s["non_gaussian_ratio_peak"] < 2.3
        assert s["lj_eps_AB"] == 0.5 and s["seed"] == 3
        assert os.path.exists(os.path.join(run_dir, "dynamic_structure.npz"))
        with open(os.path.join(run_dir, "fluctuation_summary.json")) as handle:
            written = json.load(handle)
        assert written["peak_shell_index"] == 1

    def test_frozen_pattern_gives_plateau_and_ratio_one(self, tmp_path):
        rng = np.random.default_rng(8)
        run_dir, _ = self._synthetic_run(tmp_path, rng, frozen=True)
        s = analyze_run(run_dir, write=False)["summary"]
        assert abs(s["non_gaussian_ratio_peak"] - 1.0) < 1e-6
        assert abs(s["plateau_psi_peak"] - 1.0) < 1e-6
        assert np.isnan(s["tau_psi_peak"])

    def test_window_and_lag_arguments(self, tmp_path):
        rng = np.random.default_rng(9)
        run_dir, _ = self._synthetic_run(tmp_path, rng)
        data = load_mode_amplitudes(os.path.join(run_dir, "mode_amplitudes.npz"))
        r = analyze_mode_amplitudes(data, window_fraction=0.25, max_lag_fraction=0.2)
        assert r["summary"]["n_frames_window"] == 100
        assert r["arrays"]["lag_steps"][-1] == int(0.2 * 99) * 100


@needs_openmm
class TestRunIntegration:
    def test_mode_recording_leaves_trajectory_unchanged(self, tmp_path):
        from types import SimpleNamespace

        from melt.run import execute

        base = dict(
            out=str(tmp_path), seed=5, n_chains=6, chain_length=8, box_size=6.0,
            sequence="correlated", f_A=0.5, block_length=4, kappa=0.5, pi=0.99,
            bond_k=200.0, bond_r0=1.0, lj_eps_AA=1.0, lj_eps_BB=1.0, lj_eps_AB=0.5,
            lj_eps_core=1.0, lj_sigma=1.0, lj_cutoff=2.5, temperature=1.0,
            T_equilibrate=1.0, T_quench=0.7, friction=1.0, dt=0.005,
            n_steps=400, equilibration=100, snapshot_interval=200, platform="Reference",
            compute_density=True, grid_size=8, save_trajectory=True, save_density_grids=False,
        )
        plain = execute(SimpleNamespace(run_id="plain", **base))
        modes = execute(SimpleNamespace(run_id="modes", record_modes=True, mode_interval=100, mode_q_max=2.0, **base))
        a = np.load(os.path.join(plain, "trajectory.npz"))["positions"]
        b = np.load(os.path.join(modes, "trajectory.npz"))["positions"]
        np.testing.assert_allclose(a, b, atol=1e-6)
        data = load_mode_amplitudes(os.path.join(modes, "mode_amplitudes.npz"))
        assert data["steps"].tolist() == [100, 200, 300, 400]
        assert int(data["mode_interval_steps"]) == 100
        # The recorded frame at a snapshot step reproduces the direct estimator from the trajectory.
        types = np.load(os.path.join(modes, "trajectory.npz"))["types"]
        ref = direct_partial_structure_factors(b[-1].astype(np.float64), types, data["q_vectors"])
        got = partial_structure_factors_from_amplitudes(data["rho_A"][-1], data["rho_B"][-1], types.size)
        np.testing.assert_allclose(got["S_CC"], ref["S_CC"], rtol=1e-3, atol=1e-3)
        with open(os.path.join(modes, "meta.json")) as handle:
            meta = json.load(handle)
        assert meta["mode_amplitudes"]["enabled"] is True
        assert meta["mode_amplitudes"]["n_frames"] == 4
