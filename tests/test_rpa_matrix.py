from __future__ import annotations

import json

import numpy as np
import pytest

from melt.rpa import single_chain_composition_form_factor
from melt.rpa_matrix import (balanced_intrachain_matrix, balanced_melt_rpa,
                             density_form_factor, force_equivalent_tail_transform,
                             matrix_rpa)
from scripts.submission_rpa_audit import main


def test_zero_interactions_recover_independent_chain_normalizations():
    q = np.array([0, 1e-10, 0.4, 50.0])
    omega = balanced_intrachain_matrix(q, 40, 0.7, 0.99)
    result = matrix_rpa(omega, np.zeros((2, 2)))
    s0 = single_chain_composition_form_factor(q, 40, 0.7, 0.99)
    assert np.all(result["stable"])
    assert np.allclose(result["structure"], omega)
    assert np.allclose(result["S_nn"], density_form_factor(q, 40))
    assert np.allclose(result["S_psipsi"], s0)
    assert np.allclose(result["S_npsi"], 0, atol=1e-12)
    s = result["structure"]
    assert np.allclose(s[..., 0, 0] - s[..., 0, 1], s0 / 2)
    assert np.isclose(density_form_factor(0, 40), 40)
    assert np.isclose(density_form_factor(50, 40), 1)


def test_packing_cancels_only_from_symmetric_composition():
    args = dict(chain_length=40, kappa=0.6, pi=0.99, density=0.54,
                kBT=0.7, epsilon_like=1.0, epsilon_ab=0.999)
    q = np.array([0, 0.1, 0.5, 1.0])
    a = balanced_melt_rpa(q, packing_stiffness=15, **args)
    b = balanced_melt_rpa(q, packing_stiffness=30, **args)
    s0 = single_chain_composition_form_factor(q, 40, 0.6, 0.99)
    expected_inverse = 1 / s0 + 0.54 / 0.7 * (1 - 0.999) * force_equivalent_tail_transform(q) / 2
    assert np.all(a["stable"]) and np.all(b["stable"])
    assert np.allclose(a["S_psipsi"], b["S_psipsi"])
    assert np.allclose(a["S_psipsi"], 1 / expected_inverse)
    assert np.all(b["S_nn"] < a["S_nn"])
    assert np.allclose(a["S_npsi"], 0, atol=1e-12)


def test_off_symmetry_requires_schur_complement():
    omega = np.array([[3.0, 0.4], [0.4, 1.0]])
    a = matrix_rpa(omega, np.zeros((2, 2)))
    b = matrix_rpa(omega, 2 * np.ones((2, 2)))
    for result in (a, b):
        k = result["channel_inverse"]
        assert result["stable"]
        assert np.isclose(1 / result["S_psipsi"], k[1, 1] - k[0, 1]**2 / k[0, 0])
        assert not np.isclose(1 / result["S_psipsi"], k[1, 1])
        assert not np.isclose(result["S_npsi"], 0)
    assert not np.isclose(a["S_psipsi"], b["S_psipsi"])


def test_finite_negative_reference_kernel_can_still_be_stable():
    # Effective excess inverse responses can be negative at a packing peak.
    result = balanced_melt_rpa(3, 40, 0.6, 0.99, density=0.54, kBT=0.7,
                               epsilon_like=0, epsilon_ab=0, packing_stiffness=-0.1)
    assert result["stable"]
    assert np.isclose(result["S_nn"], 1 / (1 / density_form_factor(3, 40) - 0.1))


def test_unstable_or_marginal_matrix_never_returns_negative_structure_factor():
    interaction = np.array([np.zeros((2, 2)), [[-2, 0], [0, 0]], [[-1, 0], [0, 0]]])
    result = matrix_rpa(np.eye(2), interaction)
    assert result["stable"].tolist() == [True, False, False]
    assert np.all(np.isnan(result["structure"][1:]))
    assert np.all(np.isfinite(result["inverse"]))
    # Even if composition stiffness is positive, an unstable density saddle is invalid.
    result = matrix_rpa(np.eye(2), -np.ones((2, 2)))
    assert np.isnan(result["S_psipsi"])


def test_kernel_value_small_q_and_length_units():
    q = np.array([0, 1e-12, 0.3])
    w = force_equivalent_tail_transform(q)
    assert np.isclose(w[0], -13.6581128817717, rtol=1e-12)
    assert np.isclose(w[1], w[0], rtol=1e-12)
    assert np.allclose(force_equivalent_tail_transform(q / 2, sigma=2), 8 * w)


@pytest.mark.parametrize("q,n,b", [(-0.1, 40, 1), (np.nan, 40, 1),
                                   (0, 2.5, 1), (0, 0, 1), (0, 40, 0)])
def test_invalid_chain_inputs(q, n, b):
    with pytest.raises(ValueError):
        density_form_factor(q, n, b)


@pytest.mark.parametrize("omega,interaction", [
    (np.zeros((2, 2)), np.eye(2)), (np.eye(3), np.eye(2)),
    (np.eye(2), [[1, 2], [0, 1]]), (np.eye(2), [[np.nan, 0], [0, 1]])])
def test_invalid_matrix_inputs(omega, interaction):
    with pytest.raises(ValueError):
        matrix_rpa(omega, interaction)


def test_invalid_physical_inputs():
    args = dict(chain_length=40, kappa=0.6, pi=0.99, density=0.54,
                kBT=0.7, epsilon_like=1.0, epsilon_ab=0.999, packing_stiffness=15)
    for key, value in (("kappa", -0.1), ("pi", 1.1), ("density", -1),
                       ("kBT", 0), ("packing_stiffness", np.inf)):
        with pytest.raises(ValueError):
            balanced_melt_rpa(0.1, **{**args, key: value})
    with pytest.raises(ValueError):
        force_equivalent_tail_transform(0, cutoff=1)


def test_original_submission_audit(tmp_path):
    main(["--output", str(tmp_path)])
    audit = json.loads((tmp_path / "audit.json").read_text())
    assert len(audit["modes"]) == 8
    assert np.isclose(audit["chi_bare_zero"], 9.49925635653972)
    assert audit["density_prediction"] is None
    assert audit["reference_packing_stiffness"] is None
    for mode in audit["modes"]:
        assert mode["total_bare_composition_stiffness"] < 0
        assert mode["S_psipsi_if_density_stable"] is None
    assert (tmp_path / "rpa_modes.csv").exists()
