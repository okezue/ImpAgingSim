"""Reproduce the original submission's corrected RPA diagnostics, without MD.

Run ``python -m scripts.submission_rpa_audit --output analysis/submission_rpa``.
An optional ``--source-data`` points at the submission's source_data directory
and checks the unchanged Fig. 4 descriptive regressions (legacy Fig6 CSV names).  No reference packing
stiffness is fitted or assumed, so no density structure factor is reported.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import linregress

from melt.rpa import single_chain_composition_form_factor
from melt.rpa_matrix import density_form_factor, force_equivalent_tail_transform


def _read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def _regression(points, predictor):
    x = np.asarray([1 / p[predictor] for p in points])
    y = np.asarray([1 / p["C"] for p in points])
    fit = linregress(x, y)
    return {"n_conditions": len(points), "slope": float(fit.slope),
            "intercept": float(fit.intercept), "r_squared": float(fit.rvalue**2)}


def check_source_data(source_data):
    """Recompute form factors and regressions from the supplied seed-level CSV."""
    source_data = Path(source_data)
    seeds = _read_csv(source_data / "Fig6_predictor_seed_values.csv")
    groups = {}
    errors = []
    for row in seeds:
        pi, kappa, q = (float(row[key]) for key in ("pi", "kappa", "k_star"))
        c = float(row["C_archive_N"])
        if not np.isfinite(c) or c <= 0:
            raise ValueError("source C_archive_N must be finite and positive")
        pk = float(single_chain_composition_form_factor(q, 40, kappa, pi, b=1.0))
        p0 = float(single_chain_composition_form_factor(0, 40, kappa, pi, b=1.0))
        errors.extend([abs(pk - float(row["P_N_kstar_seed"])),
                       abs(p0 - float(row["P_N_q0_seed"]))])
        groups.setdefault((pi, kappa), []).append({"C": c, "Pk": pk, "P0": p0})
    if not groups:
        raise ValueError("source seed table is empty")
    primary, inactive = [], []
    for (pi, kappa), values in groups.items():
        if pi > 0.5 and kappa > 0:
            primary.append({key: float(np.mean([v[key] for v in values]))
                            for key in ("C", "Pk", "P0")})
        else:
            inactive.extend(values)
    if len(primary) < 2 or not inactive:
        raise ValueError("source table requires active conditions and an inactive baseline")
    baseline = {key: float(np.mean([v[key] for v in inactive])) for key in ("C", "Pk", "P0")}
    fits = {name: {pred: _regression(points, pred) for pred in ("Pk", "P0")}
            for name, points in (("primary", primary), ("collapsed_boundary", primary + [baseline]))}
    max_error = float(max(errors))
    if max_error > 1e-10:
        raise ValueError(f"source form factors disagree with N=40, b=sigma=1: {max_error}")
    archived = _read_csv(source_data / "Fig6_predictor_regression_audit.csv")
    differences = []
    for row in archived:
        dataset = "collapsed_boundary" if row["dataset"].startswith("collapsed") else "primary"
        pred = "P0" if row["predictor"] == "q=0" else "Pk"
        fit = fits[dataset][pred]
        differences.extend(abs(fit[key] - float(row[key]))
                           for key in ("slope", "intercept", "r_squared"))
        if fit["n_conditions"] != int(row["n_group_means"]):
            raise ValueError("source regression condition counts do not match")
    max_fit_error = float(max(differences)) if differences else None
    if max_fit_error is None or max_fit_error > 1e-10:
        raise ValueError("recomputed regressions do not match the archived audit")
    return {"n_seed_rows": len(seeds), "n_condition_rows": len(groups),
            "max_form_factor_absolute_error": max_error,
            "max_archived_regression_absolute_error": max_fit_error,
            "fits": fits,
            "interpretation": "Response-selected descriptive fits; no chi or packing calibration."}


def build_audit():
    """Original baseline only: N=40, 144 chains, L=22 sigma, kBT=0.7."""
    n, box, rho, kbt = 40, 22.0, 5760 / 22**3, 0.7
    epsilon_like, epsilon_ab, pi = 1.0, 0.1, 0.99
    shell_norm2 = np.arange(4)
    q = 2 * np.pi / box * np.sqrt(shell_norm2)
    w = force_equivalent_tail_transform(q)
    chi = -rho * (epsilon_like - epsilon_ab) * w / kbt
    density = density_form_factor(q, n)
    rows = []
    for kappa in (0.0, 1.0):
        s0 = single_chain_composition_form_factor(q, n, kappa, pi, b=1.0)
        for i in range(len(q)):
            ideal = float(1 / s0[i])
            attraction = float(-chi[i] / 2)
            inverse = ideal + attraction
            rows.append({"kappa": kappa, "pi": pi, "shell_norm2": int(shell_norm2[i]),
                         "q_sigma": float(q[i]), "w_tilde_sigma3": float(w[i]),
                         "chi_bare": float(chi[i]), "D": float(density[i]),
                         "S0": float(s0[i]), "ideal_composition_stiffness": ideal,
                         "bare_attraction_composition_stiffness": attraction,
                         "total_bare_composition_stiffness": inverse,
                         "composition_channel_stable": inverse > 0,
                         # None, not a negative structure factor, outside the stable channel.
                         "S_psipsi_if_density_stable": float(1 / inverse) if inverse > 0 else None,
                         "mode": "formal_q0_limit" if i == 0 else "accessible_exact_shell"})
    summary = {"parameters": {"chain_length": n, "n_chains": 144, "L_sigma": box,
                              "density_sigma_minus3": rho, "kBT": kbt,
                              "epsilon_like": epsilon_like, "epsilon_ab": epsilon_ab,
                              "b_over_sigma": 1.0, "pi": pi},
               "w_tilde_zero_sigma3": float(w[0]), "chi_bare_zero": float(chi[0]),
               "reference_packing_stiffness": None,
               "density_prediction": None,
               "interpretation": [
                   "The bare attraction closure has a negative composition inverse stiffness "
                   "at every listed mode and cannot supply a physical homogeneous RPA spectrum.",
                   "This does not establish instability of the interacting MD state: the bare "
                   "closure omits correlation renormalization and is uncalibrated.",
                   "A common packing stiffness affects density but cancels from composition "
                   "only at exact ensemble A/B symmetry; no dilute limit is taken.",
                   "q=0 is a formal reference limit, not an accessible canonical mode.",
                   "No molecular-dynamics rerun is required for this audit."],
               "modes": rows}
    return summary, rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="directory for JSON and CSV outputs")
    parser.add_argument("--source-data", type=Path, help="original Supplementary Data source_data directory")
    args = parser.parse_args(argv)
    summary, rows = build_audit()
    if args.source_data is not None:
        summary["source_data_audit"] = check_source_data(args.source_data)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "rpa_modes.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (args.output / "audit.json").open("w") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(f"w_tilde(0) = {summary['w_tilde_zero_sigma3']:.12f} sigma^3")
    print(f"chi_bare(0) = {summary['chi_bare_zero']:.12f}")
    print(f"Wrote {args.output / 'rpa_modes.csv'} and {args.output / 'audit.json'}")


if __name__ == "__main__":
    main()
