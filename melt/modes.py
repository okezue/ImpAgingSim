"""Exact periodic-box Fourier amplitudes of the A and B bead densities, recorded in time.

The dynamic structure factor needs ``rho_alpha(q, t)`` at every recorded frame, not just the
final equal-time spectra kept by :mod:`melt.direct_structure`.  Because every box mode is
``q = 2*pi*h/L`` with integer ``h``, the phase factorizes per axis and one small complex
matrix product per frame gives every mode exactly, with no gridding or assignment window.
"""
from __future__ import annotations

import os
import tempfile

import numpy as np

from .direct_structure import ReciprocalModes, reciprocal_modes


class SeparableModeEvaluator:
    """Evaluate ``rho_A(q)`` and ``rho_B(q)`` for every mode in ``modes`` from bead positions."""

    def __init__(self, box_size: float, modes: ReciprocalModes) -> None:
        self.box_size = float(box_size)
        self.modes = modes
        h = np.asarray(modes.h, dtype=np.int64)
        self.h_max = int(np.max(np.abs(h))) if h.size else 0
        width = 2 * self.h_max + 1
        self._axis = np.arange(-self.h_max, self.h_max + 1, dtype=np.float64)
        shifted = h + self.h_max
        # Row index into the (hx,hy) outer product and column index into hz.
        self._row = shifted[:, 0] * width + shifted[:, 1]
        self._col = shifted[:, 2]

    def _phase_tables(self, positions: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        pos = np.asarray(positions, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"positions must have shape (N,3), got {pos.shape}")
        scale = 2.0 * np.pi / self.box_size
        ex = np.exp(1j * scale * np.outer(pos[:, 0], self._axis))
        ey = np.exp(1j * scale * np.outer(pos[:, 1], self._axis))
        ez = np.exp(1j * scale * np.outer(pos[:, 2], self._axis))
        return ex, ey, ez

    def _sum_modes(self, ex: np.ndarray, ey: np.ndarray, ez: np.ndarray) -> np.ndarray:
        if ex.shape[0] == 0:
            return np.zeros(self.modes.h.shape[0], dtype=np.complex128)
        width = ex.shape[1]
        pair = (ex[:, :, None] * ey[:, None, :]).reshape(ex.shape[0], width * width)
        cube = pair.T @ ez
        return cube[self._row, self._col]

    def evaluate(self, positions: np.ndarray, types: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        bead_types = np.asarray(types)
        pos = np.asarray(positions, dtype=np.float64)
        if bead_types.shape != (pos.shape[0],):
            raise ValueError(f"types must have shape ({pos.shape[0]},), got {bead_types.shape}")
        if not np.all(np.isin(bead_types, (0, 1))):
            raise ValueError("types must contain only 0 (B) and 1 (A)")
        ex, ey, ez = self._phase_tables(pos)
        is_A = bead_types == 1
        rho_A = self._sum_modes(ex[is_A], ey[is_A], ez[is_A])
        rho_B = self._sum_modes(ex[~is_A], ey[~is_A], ez[~is_A])
        return rho_A, rho_B


def partial_structure_factors_from_amplitudes(
    rho_A: np.ndarray, rho_B: np.ndarray, n_total: int
) -> dict[str, np.ndarray]:
    """Same normalization as :func:`melt.direct_structure.direct_partial_structure_factors`."""
    a = np.asarray(rho_A, dtype=np.complex128)
    b = np.asarray(rho_B, dtype=np.complex128)
    n = float(n_total)
    S_AA = np.abs(a) ** 2 / n
    S_BB = np.abs(b) ** 2 / n
    S_AB = np.real(a * np.conj(b)) / n
    S_CC = S_AA + S_BB - 2.0 * S_AB
    return {"S_AA": S_AA, "S_BB": S_BB, "S_AB": S_AB, "S_CC": S_CC, "S_psi_psi_over_2": 0.5 * S_CC}


class ModeAmplitudeRecorder:
    """Accumulate ``rho_A(q,t)`` and ``rho_B(q,t)`` frames and write ``mode_amplitudes.npz``."""

    def __init__(self, box_size: float, q_max: float, dtype: str = "complex64",
                 allow_type_changes: bool = False) -> None:
        self.box_size = float(box_size)
        self.q_max = float(q_max)
        self.modes = reciprocal_modes(box_size, q_max)
        self.evaluator = SeparableModeEvaluator(box_size, self.modes)
        if dtype not in ("complex64", "complex128"):
            raise ValueError(f"dtype must be complex64 or complex128, got {dtype}")
        self.dtype = np.dtype(dtype)
        self.allow_type_changes = bool(allow_type_changes)
        self.steps: list[int] = []
        self._rho_A: list[np.ndarray] = []
        self._rho_B: list[np.ndarray] = []
        self._n_A: list[int] = []
        self._n_B: list[int] = []
        self.n_A: int | None = None
        self.n_B: int | None = None

    def observe(self, step: int, positions: np.ndarray, types: np.ndarray) -> None:
        step_i = int(step)
        if self.steps and step_i <= self.steps[-1]:
            raise ValueError(f"mode frames must be recorded at strictly increasing steps, got {step_i}")
        rho_A, rho_B = self.evaluator.evaluate(positions, types)
        bead_types = np.asarray(types)
        n_A = int(np.count_nonzero(bead_types == 1))
        n_B = int(bead_types.size - n_A)
        if self.n_A is None:
            self.n_A, self.n_B = n_A, n_B
        elif (n_A, n_B) != (self.n_A, self.n_B) and not self.allow_type_changes:
            raise ValueError("bead type counts changed between frames")
        self.steps.append(step_i)
        self._n_A.append(n_A)
        self._n_B.append(n_B)
        self._rho_A.append(rho_A.astype(self.dtype))
        self._rho_B.append(rho_B.astype(self.dtype))

    @property
    def n_frames(self) -> int:
        return len(self.steps)

    def payload(self, **extra: object) -> dict[str, np.ndarray]:
        if not self.steps:
            raise RuntimeError("no mode frames were recorded")
        payload: dict[str, np.ndarray] = {
            "steps": np.asarray(self.steps, dtype=np.int64),
            "rho_A": np.asarray(self._rho_A, dtype=self.dtype),
            "rho_B": np.asarray(self._rho_B, dtype=self.dtype),
            "h": self.modes.h,
            "q_vectors": self.modes.q_vectors,
            "q": self.modes.q,
            "n2": self.modes.n2,
            "shell_n2": self.modes.shell_n2,
            "shell_q": self.modes.shell_q,
            "shell_index": self.modes.shell_index,
            "shell_degeneracy": self.modes.shell_degeneracy,
            "n_A": np.asarray(int(self.n_A)),
            "n_B": np.asarray(int(self.n_B)),
            "n_A_t": np.asarray(self._n_A, dtype=np.int64),
            "n_B_t": np.asarray(self._n_B, dtype=np.int64),
            "n_total": np.asarray(int(self.n_A + self.n_B)),
            "box_size_nm": np.asarray(self.box_size),
            "q_max_inverse_nm": np.asarray(self.q_max),
            "q_units": np.asarray("nm^-1"),
            "amplitude_definition": np.asarray("rho_alpha(q)=sum_{j in alpha} exp(i q.r_j)"),
            "normalization": np.asarray("S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total"),
            "composition_channel": np.asarray("S_psi_psi^(N)=S_CC=S_AA+S_BB-2*S_AB"),
        }
        for key, value in extra.items():
            if key in payload:
                raise ValueError(f"extra key {key!r} collides with a payload field")
            payload[key] = np.asarray(value)
        return payload

    def save(self, path: str, **extra: object) -> str:
        payload = self.payload(**extra)
        output = os.path.abspath(path)
        os.makedirs(os.path.dirname(output), exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            prefix=".mode_amplitudes_", suffix=".npz", dir=os.path.dirname(output), delete=False
        )
        temporary = handle.name
        handle.close()
        try:
            np.savez_compressed(temporary, **payload)
            os.replace(temporary, output)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return output


def load_mode_amplitudes(path: str) -> dict[str, np.ndarray]:
    """Load a ``mode_amplitudes.npz`` file into a plain dictionary of arrays."""
    with np.load(path, allow_pickle=False) as data:
        out = {key: data[key] for key in data.files}
    required = ("steps", "rho_A", "rho_B", "h", "q", "shell_index", "shell_q", "n_total", "box_size_nm")
    missing = [key for key in required if key not in out]
    if missing:
        raise ValueError(f"mode amplitude file {path} lacks fields {missing}")
    steps = out["steps"]
    if steps.ndim != 1 or out["rho_A"].shape != (steps.size, out["q"].size) or out["rho_B"].shape != out["rho_A"].shape:
        raise ValueError(f"mode amplitude file {path} has inconsistent array shapes")
    return out
