from __future__ import annotations

import numpy as np


def _unwrap(out: np.ndarray):
    """Return a Python/numpy scalar for 0-d results, the array otherwise."""
    return out[()] if out.ndim == 0 else out


def sequence_covariance(ell, kappa: float, pi: float):
    """Normalized connected covariance ``Gamma(ell)`` of the bead labels.

    ``Gamma(0)=1`` and ``Gamma(ell)=kappa**2*(2*pi-1)**ell`` for ``ell>=1``, matching
    ``melt.sequences.generate_correlated``.  Vectorized over integer ``ell>=0``.
    """
    ell_arr = np.asarray(ell)
    if not np.issubdtype(ell_arr.dtype, np.integer):
        if not np.all(np.equal(np.mod(ell_arr, 1), 0)):
            raise ValueError("ell must be integer valued")
        ell_arr = ell_arr.astype(np.int64)
    if np.any(ell_arr < 0):
        raise ValueError("ell must be non-negative")
    lam = 2.0 * float(pi) - 1.0
    # Integer exponent keeps lam**ell well defined for lam<0 (pi<0.5).
    out = np.where(ell_arr == 0, 1.0, float(kappa) ** 2 * np.power(lam, ell_arr))
    return _unwrap(np.asarray(out, dtype=np.float64))


def single_chain_composition_form_factor(q, chain_length: int, kappa: float, pi: float,
                                         f_A: float = 0.5, b: float = 1.0):
    """Per-bead composition form factor ``S0(q)`` of one ideal Gaussian chain.

    ``S0(q)=4f_A(1-f_A)[1+2*sum_{ell=1}^{N-1}(1-ell/N)*Gamma(ell)*exp(-q^2 b^2 ell/6)]``
    with ``N=chain_length`` beads and statistical segment length ``b``.  This is the
    chi=0 composition channel ``S_psipsi^(N)`` in the simulation's ``|rho_A-rho_B|^2/N_total``
    normalization, so it tends to ``4f_A(1-f_A)`` at large ``q``.  Vectorized over ``q``.
    """
    N = int(chain_length)
    if N < 1:
        raise ValueError(f"chain_length must be positive, got {chain_length}")
    q_arr = np.asarray(q, dtype=np.float64)
    ell = np.arange(1, N, dtype=np.int64)
    weights = (1.0 - ell / N) * sequence_covariance(ell, kappa, pi)
    decay = np.exp(-(q_arr[..., None] ** 2) * float(b) ** 2 * ell / 6.0)
    var = 4.0 * float(f_A) * (1.0 - float(f_A))
    out = var * (1.0 + 2.0 * np.sum(weights * decay, axis=-1))
    return _unwrap(np.asarray(out, dtype=np.float64))


def rpa_structure_factor(q, chi, chain_length: int, kappa: float, pi: float,
                         f_A: float = 0.5, b: float = 1.0):
    """Incompressible RPA composition structure factor ``S(q)=S0/(1-(chi/2)*S0)``.

    Returns NaN where the denominator is ``<=0`` (past the spinodal).  Vectorized over
    ``q`` and broadcastable ``chi``.
    """
    # With the A-indicator field the incompressible RPA is S_AA^-1 = S_AA0^-1 - 2*chi.
    # psi = 2*theta_A - 1 gives S_psipsi = 4*S_AA, hence S_psipsi^-1 = S0^-1 - chi/2.
    # Anchor: symmetric homopolymer blend has S0(0) -> N, so chi_s*N = 2 (textbook).
    S0 = np.asarray(single_chain_composition_form_factor(q, chain_length, kappa, pi, f_A, b),
                    dtype=np.float64)
    chi_arr = np.asarray(chi, dtype=np.float64)
    denom = 1.0 - 0.5 * chi_arr * S0
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(denom > 0.0, S0 / denom, np.nan)
    return _unwrap(np.asarray(out, dtype=np.float64))


def rpa_spinodal_chi(chain_length: int, kappa: float, pi: float, f_A: float = 0.5,
                     b: float = 1.0, q=None) -> tuple[float, float]:
    """Spinodal ``chi_s = 2/max_q S0(q)`` and the wavevector ``q_star`` of the maximum.

    ``S0`` is monotone decreasing in ``q`` for this sequence model, so without ``q`` the
    maximum is taken at ``q=0``.  Passing an array ``q`` restricts the search to those
    wavevectors (e.g. the box-accessible ``q>=2*pi/L``).
    """
    if q is None:
        s0_max = float(single_chain_composition_form_factor(0.0, chain_length, kappa, pi, f_A, b))
        q_star = 0.0
    else:
        q_arr = np.atleast_1d(np.asarray(q, dtype=np.float64)).ravel()
        if q_arr.size == 0:
            raise ValueError("q must contain at least one wavevector")
        s0_all = np.atleast_1d(
            single_chain_composition_form_factor(q_arr, chain_length, kappa, pi, f_A, b))
        i = int(np.nanargmax(s0_all))
        s0_max = float(s0_all[i])
        q_star = float(q_arr[i])
    chi_s = 2.0 / s0_max if s0_max > 0.0 else float("inf")
    return float(chi_s), q_star
