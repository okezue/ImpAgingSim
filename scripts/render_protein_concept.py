#!/usr/bin/env python3
"""Draw exact-count schematic copolymers and an illustrative contact topology.

    python scripts/render_protein_concept.py

Requires NumPy, SciPy, and Inkscape. Both sequence examples have 40 beads,
28 A and 12 B. Each of the three associated chains also has 40 beads: 14 A,
12 B, then 14 A. Geometry is hand-arranged, not a sampled conformation or an
equilibrium morphology. Purple dashed segments denote illustrative interchain
B contacts; solid dark segments are backbone bonds. Nothing is a membrane.

The final SVG outlines its labels for reliable display while retaining editable
vector geometry. Edit this source to change labels. PNG width is 2160 pixels.
"""

from __future__ import annotations

import html
from pathlib import Path
import subprocess

import numpy as np
from scipy.interpolate import make_interp_spline


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/figures"
A = "#007F86"
B = "#D97943"
INK = "#243345"
MUTED = "#657385"
PURPLE = "#7556A5"


def arclength_sample(points: np.ndarray, count: int) -> np.ndarray:
    length = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    t = np.linspace(0, length[-1], count)
    return np.column_stack([np.interp(t, length, points[:, d]) for d in (0, 1)])


def arm(anchors: list[tuple[float, float]]) -> np.ndarray:
    """Fifteen points: B endpoint followed by the fourteen A beads of one tail."""
    anchors = np.asarray(anchors, dtype=float)
    t = np.r_[0, np.cumsum(np.linalg.norm(np.diff(anchors, axis=0), axis=1))]
    t /= t[-1]
    raw = make_interp_spline(t, anchors, k=3)(np.linspace(0, 1, 1800))
    return arclength_sample(raw, 15)


def path(points: np.ndarray, stroke: str, width: float, extra: str = "") -> str:
    d = "M" + " L".join(f"{x:.3f},{y:.3f}" for x, y in points)
    return f'<path d="{d}" fill="none" stroke="{stroke}" stroke-width="{width}" stroke-linejoin="round" stroke-linecap="round" {extra}/>'


def circle(x: float, y: float, radius: float, color: str) -> str:
    return f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{radius}" fill="{color}" stroke="white" stroke-width="0.8"/>'


def label(x: float, y: float, content: str, size: float = 21, weight: int = 400, color: str = INK) -> str:
    return (f'<text x="{x}" y="{y}" font-family="DejaVu Sans" font-size="{size}" '
            f'font-weight="{weight}" fill="{color}">{html.escape(content)}</text>')


def chain(points: np.ndarray, sequence: np.ndarray, radius: float = 6.3) -> str:
    assert len(points) == len(sequence) == 40
    assert np.count_nonzero(sequence == 0) == 12
    assert np.count_nonzero(sequence == 1) == 28
    bond = path(points, INK, 1.7, 'opacity="0.77"')
    beads = "".join(circle(x, y, radius, A if s else B) for (x, y), s in zip(points, sequence))
    return bond + beads


def check_geometry(chains: list[np.ndarray], radius: float) -> None:
    """Verify bead visibility and prevent accidental intersections of backbones."""
    points = np.concatenate(chains)
    distances = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    np.fill_diagonal(distances, np.inf)
    pair = np.unravel_index(np.argmin(distances), distances.shape)
    assert distances.min() > 2 * radius, f"Overlapping beads: {distances.min():.3f}, {pair}, {points[list(pair)]}"

    def cross(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        d, e = b - a, c - a
        return float(d[0] * e[1] - d[1] * e[0])

    segments = [(ci, i, p[i], p[i + 1]) for ci, p in enumerate(chains) for i in range(len(p) - 1)]
    for n, (ci, i, a, b) in enumerate(segments):
        for cj, j, c, d in segments[n + 1:]:
            if ci == cj and abs(i - j) <= 1:
                continue
            intersects = cross(a, b, c) * cross(a, b, d) < 0 and cross(c, d, a) * cross(c, d, b) < 0
            assert not intersects, f"Crossing bonds: {(ci, i), (cj, j)}"


def main() -> None:
    radius = 6.3
    shapes = ['<rect width="1200" height="500" fill="white"/>']
    shapes += [label(30, 40, "A", 25, 700), label(581, 40, "B", 25, 700)]
    # The two compositions are exactly the same; only the ordering differs.
    t = np.linspace(0, 1, 2000)
    x = 53 + 435 * t
    wave = 33 * np.sin(4 * np.pi * t)
    dispersed_points = arclength_sample(np.column_stack([x, 165 + wave]), 40)
    blocky_points = arclength_sample(np.column_stack([x, 365 + wave]), 40)
    dispersed = np.ones(40, dtype=int)
    dispersed[[1, 4, 7, 11, 14, 17, 21, 24, 27, 31, 34, 37]] = 0
    blocky = np.r_[np.ones(14, int), np.zeros(12, int), np.ones(14, int)]
    check_geometry([dispersed_points, blocky_points], radius)
    shapes += [label(53, 100, "Dispersed"), label(53, 300, "Blocky")]
    shapes += [chain(dispersed_points, dispersed, radius), chain(blocky_points, blocky, radius)]
    # The only key labels are the species; every bead is countable without text.
    shapes += [circle(358, 32, 7, A), label(373, 39, "A", 19),
               circle(425, 32, 7, B), label(440, 39, "B", 19)]

    # Shared loose B contact region: three unbranched chains, each with two A tails.
    # The geometry is deliberately illustrative and does not depict a trajectory.
    cores = []
    for offset, shift in zip((874, 891, 908), (-4, 4, 0)):
        j = np.arange(12)
        cores.append(np.column_stack([offset + 7 * np.sin(2 * np.pi * j / 11),
                                      176 + 14 * j + shift]))
    upper = [
        [(874, 172), (854, 158), (822, 165), (797, 145), (773, 119), (741, 113), (713, 128), (685, 127), (660, 103)],
        [(891, 180), (894, 155), (891, 132), (874, 105), (869, 80), (892, 66), (919, 75), (944, 88)],
        [(908, 176), (931, 161), (956, 164), (978, 143), (1001, 129), (1031, 134), (1062, 157), (1093, 151), (1118, 128)],
    ]
    lower = [
        [(874, 326), (853, 343), (824, 339), (800, 357), (789, 384), (761, 400), (729, 391), (704, 364), (677, 362)],
        [(891, 334), (902, 354), (887, 376), (889, 400), (915, 415), (940, 408), (963, 425), (958, 450)],
        [(908, 330), (931, 349), (959, 347), (985, 357), (1000, 386), (1028, 399), (1058, 389), (1085, 405), (1115, 408)],
    ]
    associated = []
    for core, top, bottom in zip(cores, upper, lower):
        np.testing.assert_allclose(core[0], top[0], atol=1e-12)
        np.testing.assert_allclose(core[-1], bottom[0], atol=1e-12)
        a_top = arm(top)[1:][::-1]
        a_bottom = arm(bottom)[1:]
        associated.append(np.vstack([a_top, core, a_bottom]))
    check_geometry(associated, radius)

    # Contact dashes are distinct from covalent bonds and add no graph branch.
    for first, second, js in ((0, 1, (1, 5, 9)), (1, 2, (0, 4, 8, 11))):
        for j in js:
            shapes.append(path(np.vstack([cores[first][j], cores[second][j]]), PURPLE, 1.35,
                               'stroke-dasharray="2.2 2.3"'))
    for points in associated:
        shapes.append(chain(points, blocky, radius))

    title = "Exact-composition sequence patterns and an illustrative interchain B-contact topology"
    description = (
        "Schematic, not simulation data. Panel A: dispersed and blocky chains both contain "
        "40 beads, 28 teal A and 12 orange B. Panel B: three continuous chains with the same "
        "composition have central B segments close together and A tails fanning outward. "
        "Solid lines are backbone bonds; short purple dashes illustrate interchain contacts."
    )
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="500" '
           f'viewBox="0 0 1200 500" role="img"><title>{html.escape(title)}</title>'
           f'<desc>{html.escape(description)}</desc>{"".join(shapes)}</svg>')
    OUT.mkdir(parents=True, exist_ok=True)
    src = OUT / "protein-pattern-concept.svg"
    outlined = OUT / "protein-pattern-concept.outlined.svg"
    src.write_text(svg)
    subprocess.run(["inkscape", str(src), "--export-text-to-path", "--export-plain-svg",
                    f"--export-filename={outlined}"], check=True, capture_output=True)
    outlined.replace(src)
    subprocess.run(["inkscape", str(src), f"--export-filename={src.with_suffix('.png')}",
                    "--export-width=2160"], check=True, capture_output=True)
    print(f"Wrote {src}; verified exact bead counts, no overlaps, and no bond crossings.")


if __name__ == "__main__":
    main()
