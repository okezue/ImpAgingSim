from __future__ import annotations
import numpy as np
from melt.rpa import (sequence_covariance, single_chain_composition_form_factor,
                      rpa_structure_factor, rpa_spinodal_chi)


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
