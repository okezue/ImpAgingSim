from __future__ import annotations

from dataclasses import dataclass
import os
import tempfile

import numpy as np


@dataclass(frozen=True)
class ReciprocalModes:
    """Reciprocal-lattice modes and their cubic-shell membership."""

    h: np.ndarray
    q_vectors: np.ndarray
    q: np.ndarray
    n2: np.ndarray
    shell_n2: np.ndarray
    shell_q: np.ndarray
    shell_index: np.ndarray
    shell_degeneracy: np.ndarray


def reciprocal_modes(box_size: float, q_max: float) -> ReciprocalModes:
    """Enumerate all nonzero periodic-box modes with ``|q| <= q_max``.

    The integer triplet ``h`` represents ``q = 2*pi*h/L``.  Modes are sorted by
    ``|h|^2`` and then lexicographically, which makes files deterministic across
    runs and platforms.  A shell consists of all modes with the same ``|h|^2``.
    """

    L = float(box_size)
    q_limit = float(q_max)
    if not np.isfinite(L) or L <= 0.0:
        raise ValueError(f"box_size must be positive and finite, got {box_size}")
    if not np.isfinite(q_limit) or q_limit <= 0.0:
        raise ValueError(f"q_max must be positive and finite, got {q_max}")

    h_max = int(np.floor(q_limit * L / (2.0 * np.pi) + 1e-12))
    if h_max < 1:
        raise ValueError(
            f"q_max={q_limit} excludes the first nonzero box mode 2*pi/L={2*np.pi/L}"
        )
    axis = np.arange(-h_max, h_max + 1, dtype=np.int32)
    hx, hy, hz = np.meshgrid(axis, axis, axis, indexing="ij")
    h = np.column_stack((hx.ravel(), hy.ravel(), hz.ravel()))
    n2 = np.einsum("ij,ij->i", h, h)
    q = (2.0 * np.pi / L) * np.sqrt(n2.astype(np.float64))
    keep = (n2 > 0) & (q <= q_limit + 32.0 * np.finfo(float).eps * q_limit)
    h = h[keep]
    n2 = n2[keep]
    q = q[keep]
    order = np.lexsort((h[:, 2], h[:, 1], h[:, 0], n2))
    h = h[order]
    n2 = n2[order]
    q = q[order]

    shell_n2, shell_index, shell_degeneracy = np.unique(
        n2, return_inverse=True, return_counts=True
    )
    q_vectors = (2.0 * np.pi / L) * h.astype(np.float64)
    shell_q = (2.0 * np.pi / L) * np.sqrt(shell_n2.astype(np.float64))
    return ReciprocalModes(
        h=h,
        q_vectors=q_vectors,
        q=q,
        n2=n2,
        shell_n2=shell_n2,
        shell_q=shell_q,
        shell_index=shell_index,
        shell_degeneracy=shell_degeneracy,
    )


def direct_partial_structure_factors(
    positions: np.ndarray,
    types: np.ndarray,
    q_vectors: np.ndarray,
    chunk_size: int = 128,
) -> dict[str, np.ndarray]:
    """Compute particle-level partial structure factors for supplied modes.

    With ``rho_alpha(q) = sum_{j in alpha} exp(i q.r_j)`` and ``N`` the total
    number of beads, the returned arrays are

    ``S_AA=|rho_A|^2/N``, ``S_BB=|rho_B|^2/N``, and
    ``S_AB=Re[rho_A rho_B*]/N``.

    ``S_CC=S_AA+S_BB-2*S_AB`` is ``S_psi_psi^(N)`` for ``psi=rho_A-rho_B``.
    ``S_psi_psi_over_2=S_CC/2`` is the label-symmetric primary channel used by the
    finite-size study.  No gridding, interpolation, mean subtraction, or
    response-dependent peak selection is used.
    """

    pos = np.asarray(positions, dtype=np.float64)
    bead_types = np.asarray(types)
    qv = np.asarray(q_vectors, dtype=np.float64)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"positions must have shape (N,3), got {pos.shape}")
    if bead_types.shape != (pos.shape[0],):
        raise ValueError(
            f"types must have shape ({pos.shape[0]},), got {bead_types.shape}"
        )
    if qv.ndim != 2 or qv.shape[1] != 3:
        raise ValueError(f"q_vectors must have shape (Q,3), got {qv.shape}")
    if int(chunk_size) < 1:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")
    if not np.all(np.isin(bead_types, (0, 1))):
        raise ValueError("types must contain only 0 (B) and 1 (A)")

    n_total = pos.shape[0]
    if n_total == 0:
        raise ValueError("at least one bead is required")
    pos_A = pos[bead_types == 1]
    pos_B = pos[bead_types == 0]
    n_modes = qv.shape[0]
    S_AA = np.empty(n_modes, dtype=np.float64)
    S_BB = np.empty(n_modes, dtype=np.float64)
    S_AB = np.empty(n_modes, dtype=np.float64)

    for start in range(0, n_modes, int(chunk_size)):
        stop = min(start + int(chunk_size), n_modes)
        q_chunk = qv[start:stop]
        rho_A = np.exp(1j * (pos_A @ q_chunk.T)).sum(axis=0)
        rho_B = np.exp(1j * (pos_B @ q_chunk.T)).sum(axis=0)
        S_AA[start:stop] = np.abs(rho_A) ** 2 / n_total
        S_BB[start:stop] = np.abs(rho_B) ** 2 / n_total
        S_AB[start:stop] = np.real(rho_A * np.conj(rho_B)) / n_total

    S_CC = S_AA + S_BB - 2.0 * S_AB
    return {
        "S_AA": S_AA,
        "S_BB": S_BB,
        "S_AB": S_AB,
        "S_CC": S_CC,
        "S_psi_psi_over_2": 0.5 * S_CC,
    }


def shell_average(values: np.ndarray, shell_index: np.ndarray, n_shells: int) -> np.ndarray:
    """Average a mode-resolved one-dimensional array within cubic shells."""

    x = np.asarray(values, dtype=np.float64)
    idx = np.asarray(shell_index, dtype=np.int64)
    if x.ndim != 1 or idx.shape != x.shape:
        raise ValueError("values and shell_index must be one-dimensional arrays of equal length")
    if int(n_shells) < 1:
        raise ValueError("n_shells must be positive")
    counts = np.bincount(idx, minlength=int(n_shells)).astype(np.float64)
    if counts.size != int(n_shells) or np.any(counts == 0):
        raise ValueError("shell_index must populate every shell in [0,n_shells)")
    return np.bincount(idx, weights=x, minlength=int(n_shells)) / counts


def final_snapshot_steps(n_steps: int, snapshot_interval: int, n_frames: int) -> np.ndarray:
    """Return the final recorded production steps used for direct measurement."""

    steps = int(n_steps)
    interval = int(snapshot_interval)
    frames = int(n_frames)
    if steps < 1 or interval < 1 or frames < 1:
        raise ValueError("n_steps, snapshot_interval, and n_frames must be positive")
    n_snapshots = steps // interval
    if n_snapshots == 0:
        raise ValueError("n_steps must include at least one snapshot_interval")
    if frames > n_snapshots:
        raise ValueError(
            f"requested {frames} final frames but the schedule records only {n_snapshots} snapshots"
        )
    first = n_snapshots - frames + 1
    return interval * np.arange(first, n_snapshots + 1, dtype=np.int64)


def early_and_final_snapshot_steps(
    n_steps: int,
    snapshot_interval: int,
    early_frames: int,
    final_frames: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return adjacent earlier and final measurement blocks from a production schedule."""

    early_count = int(early_frames)
    final_count = int(final_frames)
    if early_count < 0:
        raise ValueError("early_frames must be nonnegative")
    final = final_snapshot_steps(n_steps, snapshot_interval, final_count)
    if early_count == 0:
        return np.asarray([], dtype=np.int64), final
    total = early_count + final_count
    all_steps = final_snapshot_steps(n_steps, snapshot_interval, total)
    return all_steps[:early_count], all_steps[early_count:]


class DirectStructureFactorRecorder:
    """Collect direct mode- and shell-resolved spectra at selected snapshots."""

    def __init__(
        self,
        box_size: float,
        q_max: float,
        steps: np.ndarray,
        chunk_size: int = 128,
        window_labels: np.ndarray | None = None,
    ) -> None:
        self.box_size = float(box_size)
        self.q_max = float(q_max)
        self.modes = reciprocal_modes(box_size, q_max)
        self.target_steps = np.asarray(steps, dtype=np.int64)
        if self.target_steps.ndim != 1 or self.target_steps.size == 0:
            raise ValueError("steps must be a nonempty one-dimensional array")
        if np.unique(self.target_steps).size != self.target_steps.size:
            raise ValueError("steps must not contain duplicates")
        if np.any(np.diff(self.target_steps) <= 0):
            raise ValueError("steps must be strictly increasing")
        if window_labels is None:
            self.window_labels = np.full(self.target_steps.shape, "final", dtype="U5")
        else:
            labels = np.asarray(window_labels)
            if labels.shape != self.target_steps.shape or not np.all(np.isin(labels, ("early", "final"))):
                raise ValueError("window_labels must align with steps and contain only early/final")
            self.window_labels = labels.astype("U5")
        self._target = set(int(x) for x in self.target_steps)
        self.chunk_size = int(chunk_size)
        if self.chunk_size < 1:
            raise ValueError("chunk_size must be positive")
        self.recorded_steps: list[int] = []
        self._mode_data = {
            name: [] for name in ("S_AA", "S_BB", "S_AB", "S_CC", "S_psi_psi_over_2")
        }
        self._shell_data = {
            name: []
            for name in (
                "S_AA_shell",
                "S_BB_shell",
                "S_AB_shell",
                "S_CC_shell",
                "S_psi_psi_over_2_shell",
            )
        }

    def observe(self, step: int, positions: np.ndarray, types: np.ndarray) -> bool:
        step_i = int(step)
        if step_i not in self._target:
            return False
        if step_i in self.recorded_steps:
            raise ValueError(f"direct structure factor already recorded for step {step_i}")
        spectra = direct_partial_structure_factors(
            positions, types, self.modes.q_vectors, chunk_size=self.chunk_size
        )
        self.recorded_steps.append(step_i)
        for name, values in spectra.items():
            self._mode_data[name].append(values)
            self._shell_data[f"{name}_shell"].append(
                shell_average(values, self.modes.shell_index, len(self.modes.shell_n2))
            )
        return True

    def save(self, path: str) -> str:
        if set(self.recorded_steps) != self._target:
            missing = sorted(self._target.difference(self.recorded_steps))
            raise RuntimeError(f"missing requested direct structure-factor steps: {missing}")
        order = np.argsort(np.asarray(self.recorded_steps))
        payload: dict[str, np.ndarray] = {
            "steps": np.asarray(self.recorded_steps, dtype=np.int64)[order],
            "window": self.window_labels[order],
            "h": self.modes.h,
            "q_vectors": self.modes.q_vectors,
            "q": self.modes.q,
            "n2": self.modes.n2,
            "shell_n2": self.modes.shell_n2,
            "shell_q": self.modes.shell_q,
            "shell_index": self.modes.shell_index,
            "shell_degeneracy": self.modes.shell_degeneracy,
            "box_size_nm": np.asarray(self.box_size),
            "q_max_inverse_nm": np.asarray(self.q_max),
            "q_units": np.asarray("nm^-1"),
            "normalization": np.asarray("S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total"),
            "composition_channel": np.asarray("S_psi_psi^(N)=S_CC=S_AA+S_BB-2*S_AB"),
            "primary_channel": np.asarray("S_psi_psi^(N)/2=S_CC/2"),
        }
        for name, rows in self._mode_data.items():
            payload[name] = np.asarray(rows, dtype=np.float64)[order]
        for name, rows in self._shell_data.items():
            payload[name] = np.asarray(rows, dtype=np.float64)[order]

        output = os.path.abspath(path)
        os.makedirs(os.path.dirname(output), exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            prefix=".direct_structure_factor_", suffix=".npz", dir=os.path.dirname(output), delete=False
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
