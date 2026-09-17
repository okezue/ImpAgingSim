"""Stochastic epigenetic marks: bead types that turn over and are written by local feedback.

A bead is B (marked, condensing; type 0 in this codebase) or A (unmarked; type 1).  Every
``interval_steps`` integration
steps each bead updates independently over the elapsed time ``dt_marks``:

* B -> A with rate ``k_off`` (mark turnover, the memory-loss move);
* A -> B with rate ``k_on + k_fb * H(n_B)``, where ``n_B`` is the number of B beads within
  ``r_c`` of the bead and ``H(n) = n^h / (n^h + n_half^h)`` is a saturating Hill function.
  ``k_on`` is the basal writer activity; ``k_fb`` is the reader-writer feedback gain:
  writers are recruited where readers have condensed marked beads, so an unmarked bead
  inside a B-rich region is marked faster.

Without feedback (``k_fb = 0``) the B fraction relaxes to ``k_on / (k_on + k_off)``
independently of the spatial configuration.  With feedback the local rate depends on
the structure, closing the loop between sequence and conformation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import os
import tempfile

import numpy as np
from scipy.spatial import cKDTree

TYPE_A = 1  # unmarked, open
TYPE_B = 0  # marked, condensing


@dataclass(frozen=True)
class MarkDynamics:
    k_off: float
    k_on: float
    k_fb: float = 0.0
    r_c: float = 1.5
    n_half: float = 6.0
    hill: float = 2.0
    interval_steps: int = 200

    def __post_init__(self) -> None:
        if min(self.k_off, self.k_on, self.k_fb) < 0.0:
            raise ValueError("rates must be nonnegative")
        if self.r_c <= 0.0 or self.n_half <= 0.0 or self.hill <= 0.0:
            raise ValueError("r_c, n_half and hill must be positive")
        if int(self.interval_steps) < 1:
            raise ValueError("interval_steps must be positive")

    def equilibrium_fraction_without_feedback(self) -> float:
        total = self.k_on + self.k_off
        return float(self.k_on / total) if total > 0 else float("nan")

    def as_dict(self) -> dict:
        return asdict(self)


def local_B_counts(positions: np.ndarray, types: np.ndarray, box_size: float, r_c: float) -> np.ndarray:
    """Number of B beads within ``r_c`` of every bead, periodic, excluding the bead itself."""
    pos = np.mod(np.asarray(positions, dtype=np.float64), box_size)
    is_B = np.asarray(types) == TYPE_B
    counts = np.zeros(pos.shape[0], dtype=np.int64)
    if not np.any(is_B):
        return counts
    tree_B = cKDTree(pos[is_B], boxsize=box_size)
    tree_all = cKDTree(pos, boxsize=box_size)
    neighbors = tree_all.query_ball_tree(tree_B, r_c)
    counts[:] = [len(n) for n in neighbors]
    counts[is_B] -= 1  # a B bead finds itself in the B tree
    return counts


def hill(n: np.ndarray, n_half: float, exponent: float) -> np.ndarray:
    x = np.asarray(n, dtype=np.float64) ** exponent
    return x / (x + float(n_half) ** exponent)


def step_marks(
    types: np.ndarray,
    n_B: np.ndarray,
    params: MarkDynamics,
    dt_marks: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, int, int]:
    """One stochastic update of all marks over ``dt_marks``; returns (new_types, n_on, n_off)."""
    t = np.asarray(types).astype(np.int8)
    is_B = t == TYPE_B
    p_off = 1.0 - np.exp(-params.k_off * dt_marks)
    rate_on = params.k_on + params.k_fb * hill(n_B, params.n_half, params.hill)
    p_on = 1.0 - np.exp(-rate_on * dt_marks)
    u = rng.random(t.size)
    flip_off = is_B & (u < p_off)
    flip_on = (~is_B) & (u < p_on)
    new = t.copy()
    new[flip_off] = TYPE_A
    new[flip_on] = TYPE_B
    return new, int(flip_on.sum()), int(flip_off.sum())


def mean_field_rates(params: MarkDynamics, f_B: np.ndarray, coordination: float) -> np.ndarray:
    """``df_B/dt`` for a well-mixed system with ``coordination`` neighbors per bead."""
    f = np.asarray(f_B, dtype=np.float64)
    n_B = coordination * f
    return -params.k_off * f + (1.0 - f) * (params.k_on + params.k_fb * hill(n_B, params.n_half, params.hill))


class MarkRecorder:
    """Accumulate the mark configuration at every update and write ``marks.npz``."""

    def __init__(self, params: MarkDynamics, n_beads: int) -> None:
        self.params = params
        self.n_beads = int(n_beads)
        self.steps: list[int] = []
        self._types: list[np.ndarray] = []
        self.n_on: list[int] = []
        self.n_off: list[int] = []
        self.mean_n_B: list[float] = []

    def observe(self, step: int, types: np.ndarray, n_on: int, n_off: int, n_B: np.ndarray) -> None:
        t = np.asarray(types).astype(np.int8)
        if t.shape != (self.n_beads,):
            raise ValueError(f"types must have shape ({self.n_beads},), got {t.shape}")
        self.steps.append(int(step))
        self._types.append(t.copy())
        self.n_on.append(int(n_on))
        self.n_off.append(int(n_off))
        self.mean_n_B.append(float(np.mean(n_B)))

    def payload(self, **extra: object) -> dict[str, np.ndarray]:
        if not self.steps:
            raise RuntimeError("no mark frames recorded")
        types = np.asarray(self._types, dtype=np.int8)
        payload: dict[str, np.ndarray] = {
            "steps": np.asarray(self.steps, dtype=np.int64),
            "types": types,
            "f_B": (types == TYPE_B).mean(axis=1),
            "n_on": np.asarray(self.n_on, dtype=np.int64),
            "n_off": np.asarray(self.n_off, dtype=np.int64),
            "mean_n_B": np.asarray(self.mean_n_B, dtype=np.float64),
            "type_convention": np.asarray("0 = B (marked, condensing), 1 = A (unmarked); same as the rest of the package"),
        }
        for key, value in self.params.as_dict().items():
            payload[f"param_{key}"] = np.asarray(value)
        for key, value in extra.items():
            if key in payload:
                raise ValueError(f"extra key {key!r} collides with a payload field")
            payload[key] = np.asarray(value)
        return payload

    def save(self, path: str, **extra: object) -> str:
        output = os.path.abspath(path)
        os.makedirs(os.path.dirname(output), exist_ok=True)
        handle = tempfile.NamedTemporaryFile(prefix=".marks_", suffix=".npz", dir=os.path.dirname(output), delete=False)
        temporary = handle.name
        handle.close()
        try:
            np.savez_compressed(temporary, **self.payload(**extra))
            os.replace(temporary, output)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return output


def mark_autocorrelation(types_t: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Normalized site-mark autocorrelation ``C(tau)`` averaged over sites and time origins.

    ``C(0) = 1``; for marks that forget their identity ``C`` decays to zero; a persistent
    (memorized) pattern keeps ``C`` finite.  Computed on ``m - <m>_t`` per site.
    """
    m = np.asarray(types_t, dtype=np.float64)
    T = m.shape[0]
    lag_max = T - 1 if max_lag is None else int(max_lag)
    dm = m - m.mean(axis=0, keepdims=True)
    var = float(np.mean(dm * dm))
    if var <= 0.0:
        return np.full(lag_max + 1, np.nan)
    n_fft = 1 << int(np.ceil(np.log2(2 * T)))
    spec = np.fft.rfft(dm, n=n_fft, axis=0)
    corr = np.fft.irfft(np.abs(spec) ** 2, n=n_fft, axis=0)[: lag_max + 1]
    pairs = (T - np.arange(lag_max + 1)).astype(np.float64)
    return corr.mean(axis=1) / pairs / var


def persistence_fraction(types_t: np.ndarray, value: int = TYPE_B) -> float:
    """Fraction of sites that keep ``value`` at every recorded frame."""
    m = np.asarray(types_t)
    return float(np.mean(np.all(m == value, axis=0)))
