from __future__ import annotations
import numpy as np
import pytest
from melt.rpa import (sequence_covariance, single_chain_composition_form_factor,
                      rpa_structure_factor, rpa_spinodal_chi, chi_from_eps_AB,
                      eps_AB_from_chi, segment_length_from_rg, fit_alpha_from_mixed_side)


class TestSequenceCovariance:
    def test_lag_zero_is_one(self):
        assert sequence_covariance(0, 0.5, 0.9) == 1.0
        assert sequence_covariance(np.array([0, 0]), 0.3, 0.7).tolist() == [1.0, 1.0]

    def test_lag_one_is_kappa2_lambda(self):
        kappa, pi = 0.6, 0.9
        assert np.isclose(sequence_covariance(1, kappa, pi), kappa**2 * (2 * pi - 1))

    def test_kappa_zero_kills_positive_lags(self):
        g = sequence_covariance(np.arange(0, 10), 0.0, 0.95)
        assert g[0] == 1.0
        assert np.all(g[1:] == 0.0)

    def test_vectorized_geometric_decay(self):
        kappa, pi = 1.0, 0.8
        ell = np.arange(0, 6)
        g = sequence_covariance(ell, kappa, pi)
        assert g.shape == ell.shape
        assert np.allclose(g[1:], (2 * pi - 1) ** ell[1:])

    def test_negative_lambda_integer_power(self):
        # pi<0.5 gives lambda<0; odd lags must be negative, not NaN.
        g = sequence_covariance(np.array([1, 2]), 1.0, 0.25)
        assert np.isclose(g[0], -0.5) and np.isclose(g[1], 0.25)


class TestFormFactor:
    def test_large_q_limit(self):
        for f_A in (0.5, 0.3):
            s = single_chain_composition_form_factor(50.0, 40, 1.0, 0.99, f_A=f_A, b=1.0)
            assert abs(s - 4 * f_A * (1 - f_A)) < 1e-6
        assert abs(single_chain_composition_form_factor(50.0, 40, 1.0, 0.99) - 1.0) < 1e-6

    def test_kappa_zero_flat(self):
        q = np.linspace(0.0, 5.0, 50)
        for f_A in (0.5, 0.25):
            s = single_chain_composition_form_factor(q, 40, 0.0, 0.99, f_A=f_A)
            assert s.shape == q.shape
            assert np.allclose(s, 4 * f_A * (1 - f_A))

    def test_monotone_decreasing(self):
        q = np.linspace(0.0, 10.0, 400)
        s = single_chain_composition_form_factor(q, 40, 1.0, 0.99)
        assert np.all(np.diff(s) < 0)

    def test_q_zero_matches_sum(self):
        N, kappa, pi = 40, 0.7, 0.9
        ell = np.arange(1, N)
        lam = 2 * pi - 1
        expected = 1.0 + 2.0 * np.sum((1 - ell / N) * kappa**2 * lam**ell)
        assert np.isclose(single_chain_composition_form_factor(0.0, N, kappa, pi), expected)

    def test_scalar_in_scalar_out(self):
        s = single_chain_composition_form_factor(0.3, 40, 0.5, 0.9)
        assert np.ndim(s) == 0


class TestRPA:
    def test_homopolymer_blend_anchor(self):
        N = 40
        assert np.isclose(single_chain_composition_form_factor(0.0, N, 1.0, 1.0), N)
        chi_s, q_star = rpa_spinodal_chi(N, 1.0, 1.0)
        assert np.allclose(chi_s * N, 2.0)
        assert q_star == 0.0

    def test_chi_zero_recovers_s0(self):
        q = np.linspace(0.0, 3.0, 30)
        s0 = single_chain_composition_form_factor(q, 40, 1.0, 0.99)
        assert np.allclose(rpa_structure_factor(q, 0.0, 40, 1.0, 0.99), s0)

    def test_near_and_past_spinodal(self):
        N, kappa, pi = 40, 1.0, 0.99
        chi_s, _ = rpa_spinodal_chi(N, kappa, pi)
        s0 = single_chain_composition_form_factor(0.0, N, kappa, pi)
        below = rpa_structure_factor(0.0, 0.999 * chi_s, N, kappa, pi)
        assert np.isfinite(below) and below > 100 * s0
        assert np.isnan(rpa_structure_factor(0.0, 1.01 * chi_s, N, kappa, pi))

    def test_nan_is_elementwise_over_q(self):
        N, kappa, pi = 40, 1.0, 0.99
        chi_s, _ = rpa_spinodal_chi(N, kappa, pi)
        q = np.array([0.0, 5.0])
        s = rpa_structure_factor(q, 1.5 * chi_s, N, kappa, pi)
        assert np.isnan(s[0]) and np.isfinite(s[1])

    def test_spinodal_restricted_to_q_array(self):
        N, kappa, pi = 40, 1.0, 0.99
        q = np.array([0.5, 0.3, 1.0])
        chi_s, q_star = rpa_spinodal_chi(N, kappa, pi, q=q)
        assert q_star == 0.3
        assert np.isclose(chi_s, 2.0 / single_chain_composition_form_factor(0.3, N, kappa, pi))
        assert chi_s > rpa_spinodal_chi(N, kappa, pi)[0]


class TestChiBridge:
    def test_round_trip(self):
        eps = np.array([1.0, 0.8, 0.5, 0.1])
        chi = chi_from_eps_AB(eps, 3.0, 0.7)
        assert np.allclose(eps_AB_from_chi(chi, 3.0, 0.7), eps)

    def test_zero_at_eps_like(self):
        assert chi_from_eps_AB(1.0, 3.0, 0.7) == 0.0
        assert chi_from_eps_AB(0.4, 3.0, 0.7, eps_like=0.4) == 0.0

    def test_increases_as_eps_AB_decreases(self):
        eps = np.linspace(1.0, 0.0, 11)
        chi = chi_from_eps_AB(eps, 2.0, 0.7)
        assert np.all(np.diff(chi) > 0)
        assert np.isclose(chi[-1], 2.0 / 0.7)

    def test_inverse_rejects_zero_alpha(self):
        with pytest.raises(ValueError):
            eps_AB_from_chi(0.1, 0.0, 0.7)


class TestSegmentLength:
    def test_round_trip(self):
        N, b = 40, 0.92
        rg = b * np.sqrt(N / 6.0)
        assert np.isclose(segment_length_from_rg(rg, N), b)
        rgs = np.array([1.0, 2.0])
        assert np.allclose(segment_length_from_rg(rgs, N), rgs * np.sqrt(6.0 / N))


class TestFitAlpha:
    N, kappa, pi, f_A, T_star, b = 40, 0.5, 0.99, 0.5, 0.7, 0.92
    q_peak = 2 * np.pi / 22
    alpha_true = 3.0

    def _synthetic(self, eps, seed=0, sigma=0.02):
        chi = chi_from_eps_AB(eps, self.alpha_true, self.T_star)
        S = rpa_structure_factor(self.q_peak, chi, self.N, self.kappa, self.pi, self.f_A, self.b)
        rng = np.random.default_rng(seed)
        return S * rng.lognormal(0.0, sigma, size=len(eps))

    def test_recovers_alpha(self):
        # eps_AB grid kept on the mixed side: chi_s(q_peak)~0.264 for these parameters,
        # so chi(eps_AB) must stay below that for RPA to give finite S_peak.
        eps = np.array([1.0, 0.98, 0.96, 0.95])
        S = self._synthetic(eps)
        assert np.all(np.isfinite(S))
        out = fit_alpha_from_mixed_side(eps, S, self.N, self.kappa, self.pi, self.f_A,
                                        self.T_star, self.b, self.q_peak)
        assert abs(out["alpha"] - self.alpha_true) < 0.3
        assert out["n_points"] == 4
        assert out["chi_values"].shape == (4,) and out["predicted"].shape == (4,)
        assert np.all(np.isfinite(out["predicted"]))
        assert out["residual_rms"] < 0.1
        eps_sp = out["eps_AB_spinodal"]
        assert np.isfinite(eps_sp) and 0.0 < eps_sp < 1.0
        chi_s, _ = rpa_spinodal_chi(self.N, self.kappa, self.pi, self.f_A, self.b, q=self.q_peak)
        assert np.isclose(eps_sp, eps_AB_from_chi(chi_s, out["alpha"], self.T_star))

    def test_ignores_nonfinite_points(self):
        eps = np.array([1.0, 0.98, 0.96, 0.95, 0.9])
        S = self._synthetic(eps[:4])
        S = np.concatenate([S, [np.nan]])
        out = fit_alpha_from_mixed_side(eps, S, self.N, self.kappa, self.pi, self.f_A,
                                        self.T_star, self.b, self.q_peak)
        assert out["n_points"] == 4
        assert abs(out["alpha"] - self.alpha_true) < 0.3

    def test_spec_grid_is_past_spinodal(self):
        # With alpha=3 and T*=0.7 the grid [1.0,0.9,0.8,0.7] has chi>=0.43>chi_s at
        # q_peak, so RPA yields NaN there and only one finite point survives.
        eps = np.array([1.0, 0.9, 0.8, 0.7])
        chi = chi_from_eps_AB(eps, self.alpha_true, self.T_star)
        S = rpa_structure_factor(self.q_peak, chi, self.N, self.kappa, self.pi, self.f_A, self.b)
        assert np.isfinite(S[0]) and np.all(np.isnan(S[1:]))
        with pytest.raises(ValueError):
            fit_alpha_from_mixed_side(eps, S, self.N, self.kappa, self.pi, self.f_A,
                                      self.T_star, self.b, self.q_peak)

    def test_too_few_points_raises(self):
        with pytest.raises(ValueError):
            fit_alpha_from_mixed_side([1.0, 0.9], [2.0, np.nan], self.N, self.kappa, self.pi,
                                      self.f_A, self.T_star, self.b, self.q_peak)
