"""Two-species homogeneous-reference RPA for the original submission.

All structure factors are normalized per bead.  The bare attractive-potential
closure is an uncalibrated diagnostic, not a prediction for a post-quench state.
Packing enters through an explicit effective reference stiffness; a Fourier
transform of the divergent WCA core is neither needed nor well defined here.
Existing scalar routines used by later applications are intentionally unchanged.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from .rpa import single_chain_composition_form_factor


_CHANNELS = np.array([[1.0, 1.0], [1.0, -1.0]])


def _finite(value, name):
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _chain_inputs(q, chain_length, b):
    q = _finite(q, "q")
    if np.any(q < 0):
        raise ValueError("q must be nonnegative")
    if (not np.isscalar(chain_length) or not np.isfinite(chain_length)
            or int(chain_length) != chain_length or chain_length < 1):
        raise ValueError("chain_length must be a positive integer")
    if not np.isscalar(b) or not np.isfinite(b) or b <= 0:
        raise ValueError("b must be a finite positive scalar")
    return q, int(chain_length)


def density_form_factor(q, chain_length: int, b: float = 1.0):
    """Finite Gaussian-chain density form factor, D(0)=N and D(infinity)=1."""
    q, n = _chain_inputs(q, chain_length, b)
    ell = np.arange(1, n)
    # A finite sum avoids cancellation in closed geometric-series expressions.
    return 1.0 + 2.0 * np.sum((1.0 - ell / n)
                              * np.exp(-q[..., None] ** 2 * b**2 * ell / 6), axis=-1)


def balanced_intrachain_matrix(q, chain_length: int, kappa: float, pi: float,
                               b: float = 1.0):
    """A/B intrachain covariance Ω for the balanced sequence ensemble.

    Ω=1/4 [[D+S0,D-S0],[D-S0,D+S0]].  At q=0 this is the mathematical
    long-wavelength reference, not the constrained canonical simulation mode.
    """
    q, n = _chain_inputs(q, chain_length, b)
    for value, name in ((kappa, "kappa"), (pi, "pi")):
        if not np.isscalar(value) or not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be a scalar in [0, 1]")
    d = density_form_factor(q, n, b)
    s0 = np.asarray(single_chain_composition_form_factor(q, n, kappa, pi, b=b))
    omega = np.empty(q.shape + (2, 2))
    omega[..., 0, 0] = omega[..., 1, 1] = (d + s0) / 4
    omega[..., 0, 1] = omega[..., 1, 0] = (d - s0) / 4
    return omega


def _symmetric_matrix(value, name):
    value = _finite(value, name)
    if value.ndim < 2 or value.shape[-2:] != (2, 2):
        raise ValueError(f"{name} must end in shape (2, 2)")
    if not np.allclose(value, value.swapaxes(-1, -2), rtol=1e-12, atol=1e-14):
        raise ValueError(f"{name} must be symmetric")
    return (value + value.swapaxes(-1, -2)) / 2


def matrix_rpa(omega, interaction):
    """Invert Γ=Ω^-1+V for general symmetric A/B matrices.

    ``interaction`` is the dimensionless, density-weighted effective interaction
    Hessian, including packing.  Ω must be positive definite.  All entries of
    the returned collective covariance are NaN wherever Γ is not numerically
    positive definite; unstable inverse stiffnesses remain available for audit.
    The physical channels are n=ρ_A+ρ_B and ψ=ρ_A-ρ_B, with no sqrt(2) scaling.
    Off symmetry, Sψψ^-1=Kψψ-Knψ^2/Knn, where K is ``channel_inverse``.
    """
    omega = _symmetric_matrix(omega, "omega")
    interaction = _symmetric_matrix(interaction, "interaction")
    if np.any(np.linalg.eigvalsh(omega) <= 0):
        raise ValueError("omega must be positive definite (exclude constrained zero modes)")
    inverse = np.linalg.inv(omega) + interaction
    eigenvalues = np.linalg.eigvalsh(inverse)
    tolerance = 16 * np.finfo(float).eps * np.max(np.abs(eigenvalues), axis=-1)
    stable = eigenvalues[..., 0] > tolerance
    structure = np.full(inverse.shape, np.nan)
    flat_inverse = inverse.reshape((-1, 2, 2))
    flat_structure = structure.reshape((-1, 2, 2))
    mask = np.asarray(stable).reshape(-1)
    flat_structure[mask] = np.linalg.inv(flat_inverse[mask])
    channels = _CHANNELS @ structure @ _CHANNELS.T
    channel_inverse = _CHANNELS @ inverse @ _CHANNELS.T / 4
    return {"inverse": inverse, "structure": structure, "stable": stable,
            "channel_inverse": channel_inverse,
            "S_nn": channels[..., 0, 0], "S_psipsi": channels[..., 1, 1],
            "S_npsi": channels[..., 0, 1]}


def force_equivalent_tail_transform(q, sigma: float = 1.0, cutoff: float = 2.5):
    """Spherical transform of the continuous unit-depth attractive kernel.

    ``cutoff`` is in units of sigma.  The kernel is constant below the LJ
    minimum, cutoff-shifted LJ from the minimum to cutoff, and zero outside.
    Its derivative matches the implemented attractive force on each branch.
    The constant interior is essential to avoid a spurious impulsive force.
    """
    q = _finite(q, "q")
    if np.any(q < 0):
        raise ValueError("q must be nonnegative")
    if not np.isscalar(sigma) or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be a finite positive scalar")
    if not np.isscalar(cutoff) or not np.isfinite(cutoff) or cutoff <= 2**(1 / 6):
        raise ValueError("cutoff/sigma must exceed the LJ minimum")
    rm, rc = 2**(1 / 6) * sigma, cutoff * sigma
    shift = 4 * (cutoff**-12 - cutoff**-6)

    def one(k):
        interior = quad(lambda r: r*r * (-1 - shift) * np.sinc(k*r / np.pi),
                        0, rm, epsabs=1e-11, epsrel=1e-11)[0]
        exterior = quad(lambda r: r*r * (4*((sigma/r)**12 - (sigma/r)**6) - shift)
                        * np.sinc(k*r / np.pi), rm, rc, epsabs=1e-11, epsrel=1e-11)[0]
        return 4 * np.pi * (interior + exterior)

    return np.asarray([one(k) for k in q.ravel()]).reshape(q.shape)


def balanced_melt_rpa(q, chain_length: int, kappa: float, pi: float, *,
                      density: float, kBT: float, epsilon_like: float,
                      epsilon_ab: float, packing_stiffness, b: float = 1.0,
                      sigma: float = 1.0, cutoff: float = 2.5):
    """Finite-compressibility matrix RPA with an explicitly supplied B(q).

    Γ=Ω^-1+B(q)11^T+(ρ/kBT)w_tilde(q) ε.  B is an effective reference
    packing stiffness, not the bare WCA transform; no default is assumed.
    Its Fourier-space sign is unrestricted: full Γ stability is decisive.
    At symmetry Kψψ=1/S0+ρ(ε_like-ε_ab)w_tilde/(2 kBT), whereas
    Knn=1/D+B+ρ(ε_like+ε_ab)w_tilde/(2 kBT).  The absence of B in
    Kψψ follows from mode decoupling and does not imply a dilute melt.
    """
    for value, name in ((density, "density"), (kBT, "kBT"),
                        (epsilon_like, "epsilon_like"), (epsilon_ab, "epsilon_ab")):
        if not np.isscalar(value) or not np.isfinite(value):
            raise ValueError(f"{name} must be a finite scalar")
    if density < 0 or kBT <= 0:
        raise ValueError("density must be nonnegative and kBT positive")
    packing = _finite(packing_stiffness, "packing_stiffness")
    omega = balanced_intrachain_matrix(q, chain_length, kappa, pi, b)
    w = force_equivalent_tail_transform(q, sigma, cutoff)
    epsilon = np.array([[epsilon_like, epsilon_ab], [epsilon_ab, epsilon_like]])
    interaction = packing[..., None, None] * np.ones((2, 2))
    interaction = interaction + (density / kBT * w)[..., None, None] * epsilon
    return matrix_rpa(omega, interaction)
