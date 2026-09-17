"""Condensation observables from bead positions: contact clusters of B beads and local density.

Two beads are in contact if they are within ``r_c`` under periodic boundaries.  The B
sub-network's connected components are the condensates; the summary reports how much of
the B material sits in the largest one, how many sizeable clusters exist, the B-B
coordination, and the fraction of beads in a dense environment.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .marks import TYPE_B


def contact_pairs(positions: np.ndarray, box_size: float, r_c: float) -> np.ndarray:
    pos = np.mod(np.asarray(positions, dtype=np.float64), box_size)
    tree = cKDTree(pos, boxsize=box_size)
    pairs = tree.query_pairs(r_c, output_type="ndarray")
    return pairs.reshape(-1, 2)


def cluster_labels(n: int, pairs: np.ndarray) -> np.ndarray:
    if pairs.size == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(pairs.shape[0]), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def condensation_summary(
    positions: np.ndarray,
    types: np.ndarray,
    box_size: float,
    r_c: float = 1.5,
    n_dense: int = 6,
    min_cluster: int = 5,
) -> dict[str, float]:
    pos = np.asarray(positions, dtype=np.float64)
    t = np.asarray(types)
    n = pos.shape[0]
    is_B = t == TYPE_B
    n_B = int(is_B.sum())
    pairs = contact_pairs(pos, box_size, r_c)
    same_B = is_B[pairs[:, 0]] & is_B[pairs[:, 1]]
    same_A = (~is_B[pairs[:, 0]]) & (~is_B[pairs[:, 1]])
    coordination = np.bincount(pairs.ravel(), minlength=n).astype(np.float64)
    coord_B_of_B = np.zeros(n)
    bb = pairs[same_B]
    np.add.at(coord_B_of_B, bb.ravel(), 1.0)
    out: dict[str, float] = {
        "n_B": float(n_B),
        "f_B": float(n_B / n),
        "mean_coordination": float(coordination.mean()),
        "mean_A_coordination": float(coordination[~is_B].mean()) if n_B < n else float("nan"),
        "mean_B_coordination": float(coordination[is_B].mean()) if n_B else float("nan"),
        "mean_BB_coordination": float(coord_B_of_B[is_B].mean()) if n_B else float("nan"),
        "n_BB_contacts": float(same_B.sum()),
        "n_AB_contacts": float((~same_B & ~same_A).sum()),
        "n_AA_contacts": float(same_A.sum()),
        "dense_B_fraction": float(np.mean(coord_B_of_B[is_B] >= n_dense)) if n_B else float("nan"),
        "dense_fraction": float(np.mean(coordination >= n_dense)),
    }
    if n_B:
        index_B = np.flatnonzero(is_B)
        remap = -np.ones(n, dtype=np.int64)
        remap[index_B] = np.arange(n_B)
        labels = cluster_labels(n_B, remap[bb])
        sizes = np.bincount(labels)
        out["largest_B_cluster"] = float(sizes.max())
        out["largest_B_cluster_fraction"] = float(sizes.max() / n_B)
        out["n_B_clusters"] = float(np.sum(sizes >= min_cluster))
        out["mean_B_cluster_size"] = float(np.mean(sizes[sizes >= min_cluster])) if np.any(sizes >= min_cluster) else float("nan")
        out["B_cluster_size_second_moment"] = float(np.sum(sizes.astype(np.float64) ** 2) / n_B)
    else:
        for key in ("largest_B_cluster", "largest_B_cluster_fraction", "n_B_clusters", "mean_B_cluster_size", "B_cluster_size_second_moment"):
            out[key] = float("nan")
    return out


CONDENSATION_FIELDS = (
    "step", "n_B", "f_B", "mean_coordination", "mean_A_coordination", "mean_B_coordination",
    "mean_BB_coordination", "n_BB_contacts", "n_AB_contacts", "n_AA_contacts", "dense_B_fraction",
    "dense_fraction", "largest_B_cluster", "largest_B_cluster_fraction", "n_B_clusters",
    "mean_B_cluster_size", "B_cluster_size_second_moment",
)
