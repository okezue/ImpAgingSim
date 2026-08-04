from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from statistics import multimode
from typing import Any

import numpy as np

from .direct_structure import early_and_final_snapshot_steps, reciprocal_modes, shell_average
from .integrator import KB_KJMOLK, tstar_to_kelvin
from .fixed_density_size_scan import (
    atomic_json,
    canonical_json,
    collect_provenance,
    completion_valid,
    sha256_file,
    sha256_bytes,
    validate_creation_provenance,
    validate_manifest_design_hash,
)


COMMON_Q_BIN_EDGES = np.asarray(
    [0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.50],
    dtype=np.float64,
)


def _float(value: float) -> str:
    x = float(value)
    return "" if not np.isfinite(x) else format(x, ".17g")


def _sem(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=np.float64)
    return float(np.std(x, ddof=1) / np.sqrt(x.size)) if x.size > 1 else float("nan")


def _json_float(value: float) -> float | None:
    x = float(value)
    return x if np.isfinite(x) else None


def _window_slope(steps: np.ndarray, values: np.ndarray) -> np.ndarray:
    x = np.asarray(steps, dtype=np.float64)
    y = np.asarray(values, dtype=np.float64)
    centered = x - np.mean(x)
    denominator = float(np.dot(centered, centered))
    if x.ndim != 1 or y.ndim != 2 or y.shape[0] != x.size or denominator <= 0.0:
        raise ValueError("at least two aligned frames are required for a window slope")
    return (centered @ y) / denominator


def _common_q_bin_spectra(q: np.ndarray, mode_values: np.ndarray) -> np.ndarray:
    q_values = np.asarray(q, dtype=np.float64)
    spectra = np.asarray(mode_values, dtype=np.float64)
    if spectra.ndim != 2 or spectra.shape[1] != q_values.size:
        raise ValueError("mode_values must have one column per q mode")
    columns = []
    for index, (lower, upper) in enumerate(zip(COMMON_Q_BIN_EDGES[:-1], COMMON_Q_BIN_EDGES[1:])):
        if index == len(COMMON_Q_BIN_EDGES) - 2:
            mask = (q_values >= lower) & (q_values <= upper)
        else:
            mask = (q_values >= lower) & (q_values < upper)
        if not np.any(mask):
            raise ValueError(f"common q bin [{lower},{upper}] contains no modes")
        columns.append(np.mean(spectra[:, mask], axis=1))
    return np.column_stack(columns)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _validated_spectra(
    path: Path,
    spec: dict[str, Any],
    simulation_settings: dict[str, Any],
) -> dict[str, np.ndarray]:
    required = {
        "steps",
        "window",
        "h",
        "q_vectors",
        "q",
        "n2",
        "shell_n2",
        "shell_q",
        "shell_index",
        "shell_degeneracy",
        "box_size_nm",
        "q_max_inverse_nm",
        "S_AA",
        "S_BB",
        "S_AB",
        "S_CC",
        "S_psi_psi_over_2",
        "S_AA_shell",
        "S_BB_shell",
        "S_AB_shell",
        "S_CC_shell",
        "S_psi_psi_over_2_shell",
        "normalization",
        "composition_channel",
        "primary_channel",
        "q_units",
    }
    with np.load(path, allow_pickle=False) as archive:
        missing = sorted(required.difference(archive.files))
        if missing:
            raise ValueError(f"{path} lacks arrays {missing}")
        data = {name: np.asarray(archive[name]) for name in required}
    steps = data["steps"]
    shell_q = data["shell_q"]
    shell_n2 = data["shell_n2"]
    shell_degeneracy = data["shell_degeneracy"]
    try:
        expected_early_steps, expected_final_steps = early_and_final_snapshot_steps(
            int(simulation_settings["n_steps"]),
            int(simulation_settings["snapshot_interval"]),
            int(simulation_settings["direct_early_frames"]),
            int(simulation_settings["direct_final_frames"]),
        )
        expected_steps = np.concatenate((expected_early_steps, expected_final_steps))
        q_max = float(simulation_settings["direct_q_max"])
        box_size = float(spec["box_size"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("campaign manifest lacks valid direct-spectrum settings") from exc
    if steps.dtype.kind not in "iu":
        raise ValueError(f"{path}: steps must be an integer array")
    np.testing.assert_array_equal(steps, expected_steps, err_msg=f"{path}: unexpected steps")
    expected_windows = np.asarray(
        ["early"] * len(expected_early_steps) + ["final"] * len(expected_final_steps)
    )
    np.testing.assert_array_equal(
        data["window"].astype("U5"), expected_windows, err_msg=f"{path}: window labels mismatch"
    )

    expected_modes = reciprocal_modes(box_size, q_max)
    np.testing.assert_array_equal(data["h"], expected_modes.h, err_msg=f"{path}: h mismatch")
    np.testing.assert_array_equal(data["n2"], expected_modes.n2, err_msg=f"{path}: n2 mismatch")
    np.testing.assert_array_equal(
        data["shell_n2"], expected_modes.shell_n2, err_msg=f"{path}: shell_n2 mismatch"
    )
    np.testing.assert_array_equal(
        data["shell_index"],
        expected_modes.shell_index,
        err_msg=f"{path}: shell_index mismatch",
    )
    np.testing.assert_array_equal(
        data["shell_degeneracy"],
        expected_modes.shell_degeneracy,
        err_msg=f"{path}: shell degeneracy mismatch",
    )
    float_tolerance = 64.0 * np.finfo(np.float64).eps
    np.testing.assert_allclose(
        data["q_vectors"],
        expected_modes.q_vectors,
        rtol=float_tolerance,
        atol=float_tolerance,
        err_msg=f"{path}: q_vectors mismatch",
    )
    np.testing.assert_allclose(
        data["q"],
        expected_modes.q,
        rtol=float_tolerance,
        atol=float_tolerance,
        err_msg=f"{path}: q mismatch",
    )
    np.testing.assert_allclose(
        data["shell_q"],
        expected_modes.shell_q,
        rtol=float_tolerance,
        atol=float_tolerance,
        err_msg=f"{path}: shell_q mismatch",
    )
    if data["box_size_nm"].shape != () or not np.isclose(
        float(data["box_size_nm"]), box_size, rtol=0.0, atol=float_tolerance * max(1.0, box_size)
    ):
        raise ValueError(f"{path}: box_size_nm does not match the run specification")
    if data["q_max_inverse_nm"].shape != () or not np.isclose(
        float(data["q_max_inverse_nm"]),
        q_max,
        rtol=0.0,
        atol=float_tolerance * max(1.0, q_max),
    ):
        raise ValueError(f"{path}: q_max_inverse_nm does not match the campaign manifest")

    if shell_q.ndim != 1 or shell_q.size == 0 or np.any(np.diff(shell_q) <= 0):
        raise ValueError(f"{path}: shell_q must be nonempty and strictly increasing")
    if shell_n2.shape != shell_q.shape or shell_degeneracy.shape != shell_q.shape:
        raise ValueError(f"{path}: shell metadata shapes do not agree")
    if np.any(shell_degeneracy <= 0):
        raise ValueError(f"{path}: shell degeneracies must be positive")
    expected_mode_shape = (steps.size, len(expected_modes.q))
    for name in ("S_AA", "S_BB", "S_AB", "S_CC", "S_psi_psi_over_2"):
        if data[name].shape != expected_mode_shape or not np.all(np.isfinite(data[name])):
            raise ValueError(f"{path}: {name} must be finite with shape {expected_mode_shape}")
    expected_mode_cc = data["S_AA"] + data["S_BB"] - 2.0 * data["S_AB"]
    np.testing.assert_allclose(data["S_CC"], expected_mode_cc, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(
        data["S_psi_psi_over_2"], 0.5 * expected_mode_cc, rtol=2e-13, atol=2e-13
    )

    expected_shape = (steps.size, len(expected_modes.shell_q))
    for name in (
        "S_AA_shell",
        "S_BB_shell",
        "S_AB_shell",
        "S_CC_shell",
        "S_psi_psi_over_2_shell",
    ):
        if data[name].shape != expected_shape or not np.all(np.isfinite(data[name])):
            raise ValueError(f"{path}: {name} must be finite with shape {expected_shape}")
    expected_cc = data["S_AA_shell"] + data["S_BB_shell"] - 2.0 * data["S_AB_shell"]
    np.testing.assert_allclose(data["S_CC_shell"], expected_cc, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(
        data["S_psi_psi_over_2_shell"], 0.5 * expected_cc, rtol=2e-13, atol=2e-13
    )
    for mode_name, shell_name in (
        ("S_AA", "S_AA_shell"),
        ("S_BB", "S_BB_shell"),
        ("S_AB", "S_AB_shell"),
        ("S_CC", "S_CC_shell"),
        ("S_psi_psi_over_2", "S_psi_psi_over_2_shell"),
    ):
        reconstructed = np.vstack(
            [
                shell_average(row, expected_modes.shell_index, len(expected_modes.shell_n2))
                for row in data[mode_name]
            ]
        )
        np.testing.assert_allclose(
            data[shell_name],
            reconstructed,
            rtol=2e-13,
            atol=2e-13,
            err_msg=f"{path}: {shell_name} is inconsistent with mode-resolved data",
        )
    if str(data["normalization"].item()) != "S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total":
        raise ValueError(f"{path}: unexpected normalization")
    if str(data["primary_channel"].item()) != "S_psi_psi^(N)/2=S_CC/2":
        raise ValueError(f"{path}: unexpected primary channel")
    if str(data["q_units"].item()) != "nm^-1":
        raise ValueError(f"{path}: unexpected q units")
    return data


def _validated_sequence_metadata(
    path: Path,
    spec: dict[str, Any],
    simulation_settings: dict[str, Any],
) -> dict[str, float | int]:
    with path.open() as handle:
        meta = json.load(handle)
    if not isinstance(meta, dict):
        raise ValueError(f"{path}: metadata root must be an object")
    try:
        n_beads = int(spec["n_beads"])
        sequence = np.asarray(meta["sequence"])
        melt_params = meta["melt_params"]
        run_params = meta["run_params"]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{path}: incomplete run metadata") from exc
    if sequence.shape != (n_beads,) or sequence.dtype.kind not in "biu":
        raise ValueError(f"{path}: sequence must be a binary integer vector of length {n_beads}")
    if not np.all(np.isin(sequence, (0, 1))):
        raise ValueError(f"{path}: sequence contains values outside {{0,1}}")
    if int(meta.get("n_particles", -1)) != n_beads:
        raise ValueError(f"{path}: n_particles does not match the run specification")
    if meta.get("sequence_type") != simulation_settings.get("sequence"):
        raise ValueError(f"{path}: sequence_type does not match the campaign manifest")
    if not isinstance(melt_params, dict) or not isinstance(run_params, dict):
        raise ValueError(f"{path}: melt_params and run_params must be objects")
    exact_pairs = (
        (melt_params.get("n_chains"), spec.get("n_chains"), "n_chains"),
        (melt_params.get("chain_length"), spec.get("chain_length"), "chain_length"),
        (run_params.get("n_steps"), simulation_settings.get("n_steps"), "n_steps"),
        (
            run_params.get("equilibration_steps"),
            simulation_settings.get("equilibration"),
            "equilibration_steps",
        ),
        (
            run_params.get("snapshot_interval"),
            simulation_settings.get("snapshot_interval"),
            "snapshot_interval",
        ),
        (run_params.get("seed"), spec.get("seed"), "seed"),
    )
    for actual, expected, name in exact_pairs:
        if actual != expected:
            raise ValueError(f"{path}: {name} does not match its planned value")
    try:
        actual_box = float(melt_params["box_size"])
        expected_box = float(spec["box_size"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{path}: invalid box_size") from exc
    if not np.isclose(
        actual_box,
        expected_box,
        rtol=0.0,
        atol=64.0 * np.finfo(np.float64).eps * max(1.0, abs(expected_box)),
    ):
        raise ValueError(f"{path}: box_size does not match the run specification")
    numeric_melt_settings = (
        "bond_k",
        "bond_r0",
        "lj_eps_AA",
        "lj_eps_BB",
        "lj_eps_AB",
        "lj_eps_core",
        "lj_sigma",
        "lj_cutoff",
        "friction",
        "dt",
    )
    for name in numeric_melt_settings:
        try:
            actual = float(melt_params[name])
            expected = float(simulation_settings[name])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path}: missing or invalid {name}") from exc
        if actual != expected:
            raise ValueError(f"{path}: {name} does not match the campaign manifest")
    expected_temperature_K = tstar_to_kelvin(
        simulation_settings["T_quench"], simulation_settings["lj_eps_AA"]
    )
    if float(melt_params.get("temperature", np.nan)) != expected_temperature_K:
        raise ValueError(f"{path}: melt temperature does not match T_quench")
    expected_extras = {
        "T_equilibrate_star": float(simulation_settings["T_equilibrate"]),
        "T_quench_star": float(simulation_settings["T_quench"]),
        "T_equilibrate_kelvin": tstar_to_kelvin(
            simulation_settings["T_equilibrate"], simulation_settings["lj_eps_AA"]
        ),
        "T_quench_kelvin": expected_temperature_K,
        "eps_ref_kjmol": float(simulation_settings["lj_eps_AA"]),
        "k_B_kjmolK": KB_KJMOLK,
    }
    for name, expected in expected_extras.items():
        if float(meta.get(name, np.nan)) != expected:
            raise ValueError(f"{path}: {name} does not match the campaign manifest")
    random_streams = meta.get("random_streams")
    if not isinstance(random_streams, dict) or random_streams != {
        "method": "numpy_SeedSequence_spawn",
        "root_seed": int(spec["seed"]),
        "sequence_spawn_key": [0],
        "placement_spawn_key": [1],
        "openmm_integrator_seed": int(spec["seed"]),
    }:
        raise ValueError(f"{path}: random stream metadata does not match the fixed-study scheme")
    if meta.get("sequence_parameters") != {
        "f_A": float(simulation_settings["f_A"]),
        "kappa": float(spec["kappa"]),
        "pi": float(spec["pi"]),
        "block_length": int(simulation_settings["block_length"]),
    }:
        raise ValueError(f"{path}: sequence parameters do not match the planned condition")
    if meta.get("density_analysis") != {
        "enabled": True,
        "grid_size": int(spec["grid_size"]),
    }:
        raise ValueError(f"{path}: density-analysis settings do not match the run specification")
    if meta.get("requested_openmm_platform") != simulation_settings.get("platform"):
        raise ValueError(f"{path}: requested OpenMM platform does not match the campaign manifest")
    realized_f_A = float(np.mean(sequence, dtype=np.float64))
    ensemble = meta.get("sequence_ensemble")
    if not isinstance(ensemble, dict):
        raise ValueError(f"{path}: sequence_ensemble metadata is missing")
    if simulation_settings.get("exact_global_composition") is not True:
        raise ValueError("fixed-density campaign must enable exact_global_composition")
    expected_method = "independent_per_chain_rejection_conditioned_on_exact_global_count"
    if ensemble.get("method") != expected_method:
        raise ValueError(f"{path}: unexpected sequence-ensemble method")
    max_attempts = simulation_settings.get("sequence_balance_max_attempts")
    attempts = ensemble.get("acceptance_attempts")
    if (
        not isinstance(max_attempts, int)
        or max_attempts < 1
        or not isinstance(attempts, int)
        or attempts < 1
        or attempts > max_attempts
        or ensemble.get("max_attempts") != max_attempts
    ):
        raise ValueError(f"{path}: invalid sequence-balance attempt metadata")
    target_A = n_beads // 2
    realized_A = int(np.sum(sequence, dtype=np.int64))
    if (
        n_beads % 2
        or realized_A != target_A
        or ensemble.get("target_A_count") != target_A
        or ensemble.get("realized_A_count") != realized_A
        or ensemble.get("constructed_complement_pairs") is not False
        or not np.isclose(float(ensemble.get("exact_global_f_A", np.nan)), 0.5)
        or realized_f_A != 0.5
    ):
        raise ValueError(f"{path}: sequence ensemble is not exactly globally balanced")
    return {"realized_f_A": realized_f_A, "sequence_balance_attempts": attempts}


def _validate_planned_run(
    path: Path,
    spec: dict[str, Any],
    simulation_settings: dict[str, Any],
) -> None:
    with path.open() as handle:
        planned = json.load(handle)
    if not isinstance(planned, dict):
        raise ValueError(f"{path}: planned run must be an object")
    if planned.get("run_spec") != spec:
        raise ValueError(f"{path}: run_spec does not match the campaign manifest")
    if planned.get("simulation_settings") != simulation_settings:
        raise ValueError(f"{path}: simulation_settings do not match the campaign manifest")


def aggregate_campaign(campaign_dir: Path | str, require_complete: bool = True) -> dict[str, Path]:
    """Aggregate final-frame direct spectra with deterministic ordering and formulas."""

    root = Path(campaign_dir).resolve()
    manifest_path = root / "manifest.json"
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    validate_manifest_design_hash(manifest)
    repo_root = Path(__file__).resolve().parents[1]
    analysis_provenance = collect_provenance(repo_root)
    campaign_commit = validate_creation_provenance(manifest, analysis_provenance)
    simulation_settings = manifest.get("simulation_settings")
    if not isinstance(simulation_settings, dict):
        raise ValueError("campaign manifest lacks simulation_settings")
    required_settings = {
        "sequence",
        "f_A",
        "block_length",
        "bond_k",
        "bond_r0",
        "lj_eps_AA",
        "lj_eps_BB",
        "lj_eps_AB",
        "lj_eps_core",
        "lj_sigma",
        "lj_cutoff",
        "T_equilibrate",
        "T_quench",
        "friction",
        "dt",
        "platform",
        "compute_density",
        "compute_direct_structure_factor",
        "exact_global_composition",
        "sequence_balance_max_attempts",
        "equilibration",
        "n_steps",
        "snapshot_interval",
        "direct_final_frames",
        "direct_early_frames",
        "direct_q_max",
        "save_trajectory",
    }
    missing_settings = sorted(required_settings.difference(simulation_settings))
    if missing_settings:
        raise ValueError(f"campaign manifest lacks simulation settings {missing_settings}")
    if not isinstance(simulation_settings["save_trajectory"], bool):
        raise ValueError("campaign save_trajectory setting must be boolean")
    if simulation_settings["compute_density"] is not True or simulation_settings[
        "compute_direct_structure_factor"
    ] is not True:
        raise ValueError("fixed-density campaign must enable both gridded and direct spectra")
    if float(simulation_settings["f_A"]) != 0.5:
        raise ValueError("fixed-density campaign must use f_A=0.5")
    if float(simulation_settings["lj_eps_AA"]) != float(simulation_settings["lj_eps_BB"]):
        raise ValueError("fixed-density campaign must use lj_eps_AA=lj_eps_BB")
    if (
        int(simulation_settings["direct_early_frames"]) != 5
        or int(simulation_settings["direct_final_frames"]) != 5
        or float(simulation_settings["direct_q_max"]) != 1.5
    ):
        raise ValueError("fixed-density campaign must use 5+5 direct frames and q_max=1.5 nm^-1")
    if int(simulation_settings["n_steps"]) % int(simulation_settings["snapshot_interval"]) != 0:
        raise ValueError("campaign n_steps is not divisible by snapshot_interval")
    direct_design = manifest.get("direct_structure_factor")
    if not isinstance(direct_design, dict) or direct_design.get(
        "common_q_bin_edges_inverse_nm"
    ) != COMMON_Q_BIN_EDGES.tolist():
        raise ValueError("campaign common-q-bin design does not match the registered analysis")
    plan = sorted(manifest["runs"], key=lambda row: int(row["index"]))
    run_records: list[dict[str, Any]] = []
    missing_runs: list[str] = []
    source_files: dict[str, str] = {}

    for spec in plan:
        run_id = str(spec["run_id"])
        run_dir = root / "runs" / run_id
        if not completion_valid(
            run_dir,
            spec,
            save_trajectory=simulation_settings["save_trajectory"],
            expected_git_commit=campaign_commit,
        ):
            missing_runs.append(run_id)
            continue
        direct_path = run_dir / "direct_structure_factor.npz"
        planned_path = run_dir / "planned_run.json"
        _validate_planned_run(planned_path, spec, simulation_settings)
        data = _validated_spectra(direct_path, spec, simulation_settings)
        meta_path = run_dir / "meta.json"
        sequence_metadata = _validated_sequence_metadata(meta_path, spec, simulation_settings)
        realized_f_A = float(sequence_metadata["realized_f_A"])
        sequence_balance_attempts = int(sequence_metadata["sequence_balance_attempts"])
        source_files[str(direct_path.relative_to(root))] = sha256_file(direct_path)
        source_files[str(meta_path.relative_to(root))] = sha256_file(meta_path)
        source_files[str(planned_path.relative_to(root))] = sha256_file(planned_path)
        early_mask = data["window"].astype("U5") == "early"
        final_mask = data["window"].astype("U5") == "final"
        if int(np.sum(early_mask)) != int(simulation_settings["direct_early_frames"]):
            raise ValueError(f"{direct_path}: wrong number of early frames")
        if int(np.sum(final_mask)) != int(simulation_settings["direct_final_frames"]):
            raise ValueError(f"{direct_path}: wrong number of final frames")
        shell_final_means = {
            name: np.mean(data[name][final_mask], axis=0)
            for name in (
                "S_AA_shell",
                "S_BB_shell",
                "S_AB_shell",
                "S_CC_shell",
                "S_psi_psi_over_2_shell",
            )
        }
        shell_early_means = {
            name: np.mean(data[name][early_mask], axis=0)
            for name in (
                "S_AA_shell",
                "S_BB_shell",
                "S_AB_shell",
                "S_CC_shell",
                "S_psi_psi_over_2_shell",
            )
        }
        primary = shell_final_means["S_psi_psi_over_2_shell"]
        primary_early = shell_early_means["S_psi_psi_over_2_shell"]
        primary_change = primary - primary_early
        primary_slope = _window_slope(
            data["steps"][final_mask], data["S_psi_psi_over_2_shell"][final_mask]
        )
        diagnostic = shell_final_means["S_AA_shell"] - shell_final_means["S_AB_shell"]
        common_frame_spectra = _common_q_bin_spectra(data["q"], data["S_psi_psi_over_2"])
        common_final = np.mean(common_frame_spectra[final_mask], axis=0)
        common_early = np.mean(common_frame_spectra[early_mask], axis=0)
        common_change = common_final - common_early
        common_slope = _window_slope(data["steps"][final_mask], common_frame_spectra[final_mask])
        peak_index = int(np.argmax(primary))
        diagnostic_peak_index = int(np.argmax(diagnostic))
        run_records.append(
            {
                "spec": spec,
                "n_frames": int(data["steps"].size),
                "n_early_frames": int(np.sum(early_mask)),
                "n_final_frames": int(np.sum(final_mask)),
                "early_first_step": int(data["steps"][early_mask][0]),
                "early_last_step": int(data["steps"][early_mask][-1]),
                "final_first_step": int(data["steps"][final_mask][0]),
                "final_last_step": int(data["steps"][final_mask][-1]),
                "first_step": int(data["steps"][0]),
                "last_step": int(data["steps"][-1]),
                "shell_q": data["shell_q"].astype(np.float64),
                "shell_n2": data["shell_n2"].astype(np.int64),
                "shell_degeneracy": data["shell_degeneracy"].astype(np.int64),
                "spectra": shell_final_means,
                "early_spectra": shell_early_means,
                "primary": primary,
                "primary_early": primary_early,
                "primary_change": primary_change,
                "primary_slope": primary_slope,
                "common_primary": common_final,
                "common_primary_early": common_early,
                "common_primary_change": common_change,
                "common_primary_slope": common_slope,
                "diagnostic": diagnostic,
                "peak_index": peak_index,
                "k_star": float(data["shell_q"][peak_index]),
                "peak_amplitude": float(primary[peak_index]),
                "realized_f_A": realized_f_A,
                "sequence_balance_attempts": sequence_balance_attempts,
                "diagnostic_peak_index": diagnostic_peak_index,
                "diagnostic_k_star": float(data["shell_q"][diagnostic_peak_index]),
                "diagnostic_peak_amplitude": float(diagnostic[diagnostic_peak_index]),
            }
        )

    if missing_runs and require_complete:
        raise RuntimeError(
            f"{len(missing_runs)} of {len(plan)} planned runs are missing or invalid; "
            "rerun them or use --allow-incomplete-analysis"
        )
    if not run_records:
        raise RuntimeError("no complete runs are available for analysis")

    per_seed_rows: list[dict[str, Any]] = []
    shell_seed_rows: list[dict[str, Any]] = []
    for record in run_records:
        spec = record["spec"]
        per_seed_rows.append(
            {
                "run_index": int(spec["index"]),
                "run_id": spec["run_id"],
                "n_chains": int(spec["n_chains"]),
                "chain_length": int(spec["chain_length"]),
                "n_beads": int(spec["n_beads"]),
                "box_size": _float(spec["box_size"]),
                "bead_density": _float(spec["bead_density"]),
                "kappa": _float(spec["kappa"]),
                "pi": _float(spec["pi"]),
                "seed": int(spec["seed"]),
                "realized_f_A": _float(record["realized_f_A"]),
                "sequence_balance_attempts": int(record["sequence_balance_attempts"]),
                "n_early_frames": record["n_early_frames"],
                "n_final_frames": record["n_final_frames"],
                "early_first_step": record["early_first_step"],
                "early_last_step": record["early_last_step"],
                "final_first_step": record["final_first_step"],
                "final_last_step": record["final_last_step"],
                "secondary_per_seed_max_metric": (
                    "max_exact_shell mean_final_frames[S_psi_psi^(N)(q)/2]"
                ),
                "secondary_per_seed_max_k_star": _float(record["k_star"]),
                "secondary_per_seed_max_k_star_units": "nm^-1",
                "secondary_per_seed_max_amplitude": _float(record["peak_amplitude"]),
                "secondary_per_seed_max_boundary": (
                    "q_min"
                    if record["peak_index"] == 0
                    else "q_max"
                    if record["peak_index"] == len(record["shell_q"]) - 1
                    else "none"
                ),
                "secondary_per_seed_max_early_to_final_change": _float(
                    record["primary_change"][record["peak_index"]]
                ),
                "secondary_per_seed_max_final_window_slope_per_step": _float(
                    record["primary_slope"][record["peak_index"]]
                ),
                "diagnostic_AA_minus_AB_at_secondary_per_seed_max_k_star": _float(
                    record["diagnostic"][record["peak_index"]]
                ),
                "diagnostic_AA_minus_AB_k_star": _float(record["diagnostic_k_star"]),
                "diagnostic_AA_minus_AB_peak_amplitude": _float(
                    record["diagnostic_peak_amplitude"]
                ),
            }
        )
        for shell_index in range(len(record["shell_q"])):
            shell_seed_rows.append(
                {
                    "run_index": int(spec["index"]),
                    "run_id": spec["run_id"],
                    "n_chains": int(spec["n_chains"]),
                    "kappa": _float(spec["kappa"]),
                    "pi": _float(spec["pi"]),
                    "seed": int(spec["seed"]),
                    "shell_index": shell_index,
                    "shell_n2": int(record["shell_n2"][shell_index]),
                    "degeneracy": int(record["shell_degeneracy"][shell_index]),
                    "q": _float(record["shell_q"][shell_index]),
                    "q_units": "nm^-1",
                    "S_AA": _float(record["spectra"]["S_AA_shell"][shell_index]),
                    "S_BB": _float(record["spectra"]["S_BB_shell"][shell_index]),
                    "S_AB": _float(record["spectra"]["S_AB_shell"][shell_index]),
                    "S_CC": _float(record["spectra"]["S_CC_shell"][shell_index]),
                    "S_psi_psi_over_2": _float(record["primary"][shell_index]),
                    "S_psi_psi_over_2_early": _float(record["primary_early"][shell_index]),
                    "S_psi_psi_over_2_early_to_final_change": _float(
                        record["primary_change"][shell_index]
                    ),
                    "S_psi_psi_over_2_final_window_slope_per_step": _float(
                        record["primary_slope"][shell_index]
                    ),
                    "diagnostic_AA_minus_AB": _float(record["diagnostic"][shell_index]),
                }
            )

    groups: dict[tuple[int, float], list[dict[str, Any]]] = {}
    per_seed_by_run = {row["run_id"]: row for row in per_seed_rows}
    for record in run_records:
        spec = record["spec"]
        groups.setdefault((int(spec["n_chains"]), float(spec["kappa"])), []).append(record)

    condition_rows: list[dict[str, Any]] = []
    shell_condition_rows: list[dict[str, Any]] = []
    common_bin_condition_rows: list[dict[str, Any]] = []
    json_conditions: list[dict[str, Any]] = []
    for key in sorted(groups):
        records = sorted(groups[key], key=lambda row: int(row["spec"]["seed"]))
        reference_q = records[0]["shell_q"]
        reference_n2 = records[0]["shell_n2"]
        reference_deg = records[0]["shell_degeneracy"]
        for record in records[1:]:
            np.testing.assert_array_equal(record["shell_n2"], reference_n2)
            np.testing.assert_allclose(record["shell_q"], reference_q, rtol=0.0, atol=0.0)
            np.testing.assert_array_equal(record["shell_degeneracy"], reference_deg)

        secondary_peak_values = np.asarray([record["peak_amplitude"] for record in records])
        secondary_k_values = np.asarray([record["k_star"] for record in records])
        primary_matrix = np.vstack([record["primary"] for record in records])
        mean_primary = np.mean(primary_matrix, axis=0)
        mean_peak_index = int(np.argmax(mean_primary))
        selected_values = primary_matrix[:, mean_peak_index]
        selected_changes = np.asarray(
            [record["primary_change"][mean_peak_index] for record in records]
        )
        selected_slopes = np.asarray(
            [record["primary_slope"][mean_peak_index] for record in records]
        )
        exact_boundary = (
            "q_min"
            if mean_peak_index == 0
            else "q_max"
            if mean_peak_index == len(reference_q) - 1
            else "none"
        )
        common_matrix = np.vstack([record["common_primary"] for record in records])
        common_mean = np.mean(common_matrix, axis=0)
        common_peak_index = int(np.argmax(common_mean))
        common_selected_values = common_matrix[:, common_peak_index]
        common_selected_changes = np.asarray(
            [record["common_primary_change"][common_peak_index] for record in records]
        )
        common_selected_slopes = np.asarray(
            [record["common_primary_slope"][common_peak_index] for record in records]
        )
        common_boundary = (
            "q_min"
            if common_peak_index == 0
            else "q_max"
            if common_peak_index == len(COMMON_Q_BIN_EDGES) - 2
            else "none"
        )
        diagnostic_peak_values = np.asarray(
            [record["diagnostic_peak_amplitude"] for record in records]
        )
        diagnostic_k_values = np.asarray([record["diagnostic_k_star"] for record in records])
        realized_f_A_values = np.asarray([record["realized_f_A"] for record in records])
        balance_attempt_values = np.asarray(
            [record["sequence_balance_attempts"] for record in records], dtype=np.float64
        )
        modes = sorted(multimode([float(value) for value in secondary_k_values]))
        modal_k = modes[0]
        spec0 = records[0]["spec"]
        for record in records:
            seed_row = per_seed_by_run[record["spec"]["run_id"]]
            seed_row.update(
                {
                    "primary_condition_selected_shell_index": mean_peak_index,
                    "primary_condition_selected_k_star": _float(reference_q[mean_peak_index]),
                    "primary_condition_selected_k_star_units": "nm^-1",
                    "primary_condition_selected_boundary": exact_boundary,
                    "primary_at_condition_selected_shell_final": _float(
                        record["primary"][mean_peak_index]
                    ),
                    "primary_at_condition_selected_shell_early": _float(
                        record["primary_early"][mean_peak_index]
                    ),
                    "primary_at_condition_selected_shell_early_to_final_change": _float(
                        record["primary_change"][mean_peak_index]
                    ),
                    "primary_at_condition_selected_shell_final_window_slope_per_step": _float(
                        record["primary_slope"][mean_peak_index]
                    ),
                }
            )
        condition = {
            "n_chains": key[0],
            "chain_length": int(spec0["chain_length"]),
            "n_beads": int(spec0["n_beads"]),
            "box_size": _float(spec0["box_size"]),
            "bead_density": _float(spec0["bead_density"]),
            "kappa": _float(key[1]),
            "pi": _float(spec0["pi"]),
            "n_seeds": len(records),
            "realized_f_A_mean": _float(np.mean(realized_f_A_values)),
            "realized_f_A_sem": _float(_sem(realized_f_A_values)),
            "sequence_balance_attempts_mean": _float(np.mean(balance_attempt_values)),
            "sequence_balance_attempts_sem": _float(_sem(balance_attempt_values)),
            "primary_selection_rule": "maximum of across-seed mean final-window exact-shell spectrum",
            "primary_candidate_shell_count": len(reference_q),
            "primary_selected_shell_index": mean_peak_index,
            "primary_k_star": _float(reference_q[mean_peak_index]),
            "primary_k_star_units": "nm^-1",
            "primary_boundary": exact_boundary,
            "primary_amplitude_seed_mean_at_selected_shell": _float(np.mean(selected_values)),
            "primary_amplitude_seed_sem_at_selected_shell": _float(_sem(selected_values)),
            "primary_early_to_final_change_seed_mean_at_selected_shell": _float(
                np.mean(selected_changes)
            ),
            "primary_early_to_final_change_seed_sem_at_selected_shell": _float(
                _sem(selected_changes)
            ),
            "primary_final_window_slope_seed_mean_per_step_at_selected_shell": _float(
                np.mean(selected_slopes)
            ),
            "primary_final_window_slope_seed_sem_per_step_at_selected_shell": _float(
                _sem(selected_slopes)
            ),
            "secondary_per_seed_max_amplitude_mean": _float(np.mean(secondary_peak_values)),
            "secondary_per_seed_max_amplitude_sem": _float(_sem(secondary_peak_values)),
            "secondary_per_seed_max_k_star_mean": _float(np.mean(secondary_k_values)),
            "secondary_per_seed_max_k_star_sem": _float(_sem(secondary_k_values)),
            "secondary_per_seed_max_k_star_mode": _float(modal_k),
            "common_bin_candidate_count": len(COMMON_Q_BIN_EDGES) - 1,
            "common_bin_selected_index": common_peak_index,
            "common_bin_selected_q_low": _float(COMMON_Q_BIN_EDGES[common_peak_index]),
            "common_bin_selected_q_high": _float(COMMON_Q_BIN_EDGES[common_peak_index + 1]),
            "common_bin_selected_boundary": common_boundary,
            "common_bin_amplitude_seed_mean": _float(np.mean(common_selected_values)),
            "common_bin_amplitude_seed_sem": _float(_sem(common_selected_values)),
            "common_bin_early_to_final_change_seed_mean": _float(
                np.mean(common_selected_changes)
            ),
            "common_bin_early_to_final_change_seed_sem": _float(_sem(common_selected_changes)),
            "common_bin_final_window_slope_seed_mean_per_step": _float(
                np.mean(common_selected_slopes)
            ),
            "common_bin_final_window_slope_seed_sem_per_step": _float(
                _sem(common_selected_slopes)
            ),
            "diagnostic_AA_minus_AB_peak_seed_mean": _float(np.mean(diagnostic_peak_values)),
            "diagnostic_AA_minus_AB_peak_seed_sem": _float(_sem(diagnostic_peak_values)),
            "diagnostic_AA_minus_AB_k_star_seed_mean": _float(np.mean(diagnostic_k_values)),
            "diagnostic_AA_minus_AB_k_star_seed_sem": _float(_sem(diagnostic_k_values)),
        }
        condition_rows.append(condition)
        json_conditions.append(
            {
                "n_chains": key[0],
                "chain_length": int(spec0["chain_length"]),
                "n_beads": int(spec0["n_beads"]),
                "box_size": float(spec0["box_size"]),
                "bead_density": float(spec0["bead_density"]),
                "kappa": key[1],
                "pi": float(spec0["pi"]),
                "n_seeds": len(records),
                "realized_f_A_mean": _json_float(np.mean(realized_f_A_values)),
                "realized_f_A_sem": _json_float(_sem(realized_f_A_values)),
                "sequence_balance_attempts_mean": _json_float(np.mean(balance_attempt_values)),
                "sequence_balance_attempts_sem": _json_float(_sem(balance_attempt_values)),
                "primary_selection_rule": condition["primary_selection_rule"],
                "primary_candidate_shell_count": len(reference_q),
                "primary_selected_shell_index": mean_peak_index,
                "primary_k_star": _json_float(reference_q[mean_peak_index]),
                "primary_k_star_units": "nm^-1",
                "primary_boundary": exact_boundary,
                "primary_amplitude_seed_mean_at_selected_shell": _json_float(
                    np.mean(selected_values)
                ),
                "primary_amplitude_seed_sem_at_selected_shell": _json_float(
                    _sem(selected_values)
                ),
                "primary_early_to_final_change_seed_mean_at_selected_shell": _json_float(
                    np.mean(selected_changes)
                ),
                "primary_early_to_final_change_seed_sem_at_selected_shell": _json_float(
                    _sem(selected_changes)
                ),
                "primary_final_window_slope_seed_mean_per_step_at_selected_shell": _json_float(
                    np.mean(selected_slopes)
                ),
                "primary_final_window_slope_seed_sem_per_step_at_selected_shell": _json_float(
                    _sem(selected_slopes)
                ),
                "secondary_per_seed_max_amplitude_mean": _json_float(
                    np.mean(secondary_peak_values)
                ),
                "secondary_per_seed_max_amplitude_sem": _json_float(
                    _sem(secondary_peak_values)
                ),
                "secondary_per_seed_max_k_star_mean": _json_float(
                    np.mean(secondary_k_values)
                ),
                "secondary_per_seed_max_k_star_sem": _json_float(_sem(secondary_k_values)),
                "secondary_per_seed_max_k_star_mode": _json_float(modal_k),
                "common_bin_candidate_count": len(COMMON_Q_BIN_EDGES) - 1,
                "common_bin_selected_index": common_peak_index,
                "common_bin_selected_q_low": _json_float(
                    COMMON_Q_BIN_EDGES[common_peak_index]
                ),
                "common_bin_selected_q_high": _json_float(
                    COMMON_Q_BIN_EDGES[common_peak_index + 1]
                ),
                "common_bin_selected_boundary": common_boundary,
                "common_bin_amplitude_seed_mean": _json_float(
                    np.mean(common_selected_values)
                ),
                "common_bin_amplitude_seed_sem": _json_float(_sem(common_selected_values)),
                "common_bin_early_to_final_change_seed_mean": _json_float(
                    np.mean(common_selected_changes)
                ),
                "common_bin_early_to_final_change_seed_sem": _json_float(
                    _sem(common_selected_changes)
                ),
                "common_bin_final_window_slope_seed_mean_per_step": _json_float(
                    np.mean(common_selected_slopes)
                ),
                "common_bin_final_window_slope_seed_sem_per_step": _json_float(
                    _sem(common_selected_slopes)
                ),
                "diagnostic_AA_minus_AB_peak_seed_mean": _json_float(
                    np.mean(diagnostic_peak_values)
                ),
                "diagnostic_AA_minus_AB_peak_seed_sem": _json_float(
                    _sem(diagnostic_peak_values)
                ),
                "diagnostic_AA_minus_AB_k_star_seed_mean": _json_float(
                    np.mean(diagnostic_k_values)
                ),
                "diagnostic_AA_minus_AB_k_star_seed_sem": _json_float(
                    _sem(diagnostic_k_values)
                ),
            }
        )

        for shell_index in range(len(reference_q)):
            row: dict[str, Any] = {
                "n_chains": key[0],
                "kappa": _float(key[1]),
                "pi": _float(spec0["pi"]),
                "n_seeds": len(records),
                "shell_index": shell_index,
                "shell_n2": int(reference_n2[shell_index]),
                "degeneracy": int(reference_deg[shell_index]),
                "q": _float(reference_q[shell_index]),
                "q_units": "nm^-1",
                "selected_as_primary": shell_index == mean_peak_index,
            }
            for output_name, array_name in (
                ("S_AA", "S_AA_shell"),
                ("S_BB", "S_BB_shell"),
                ("S_AB", "S_AB_shell"),
                ("S_CC", "S_CC_shell"),
                ("S_psi_psi_over_2", "S_psi_psi_over_2_shell"),
            ):
                values = np.asarray([record["spectra"][array_name][shell_index] for record in records])
                row[f"{output_name}_mean"] = _float(np.mean(values))
                row[f"{output_name}_sem"] = _float(_sem(values))
            diagnostic_values = np.asarray(
                [record["diagnostic"][shell_index] for record in records]
            )
            row["diagnostic_AA_minus_AB_mean"] = _float(np.mean(diagnostic_values))
            row["diagnostic_AA_minus_AB_sem"] = _float(_sem(diagnostic_values))
            early_values = np.asarray(
                [record["primary_early"][shell_index] for record in records]
            )
            change_values = np.asarray(
                [record["primary_change"][shell_index] for record in records]
            )
            slope_values = np.asarray(
                [record["primary_slope"][shell_index] for record in records]
            )
            row["S_psi_psi_over_2_early_mean"] = _float(np.mean(early_values))
            row["S_psi_psi_over_2_early_sem"] = _float(_sem(early_values))
            row["S_psi_psi_over_2_early_to_final_change_mean"] = _float(
                np.mean(change_values)
            )
            row["S_psi_psi_over_2_early_to_final_change_sem"] = _float(
                _sem(change_values)
            )
            row["S_psi_psi_over_2_final_window_slope_mean_per_step"] = _float(
                np.mean(slope_values)
            )
            row["S_psi_psi_over_2_final_window_slope_sem_per_step"] = _float(
                _sem(slope_values)
            )
            shell_condition_rows.append(row)

        for bin_index in range(len(COMMON_Q_BIN_EDGES) - 1):
            bin_values = np.asarray([record["common_primary"][bin_index] for record in records])
            bin_early = np.asarray(
                [record["common_primary_early"][bin_index] for record in records]
            )
            bin_change = np.asarray(
                [record["common_primary_change"][bin_index] for record in records]
            )
            bin_slope = np.asarray(
                [record["common_primary_slope"][bin_index] for record in records]
            )
            common_bin_condition_rows.append(
                {
                    "n_chains": key[0],
                    "kappa": _float(key[1]),
                    "pi": _float(spec0["pi"]),
                    "n_seeds": len(records),
                    "bin_index": bin_index,
                    "q_low": _float(COMMON_Q_BIN_EDGES[bin_index]),
                    "q_high": _float(COMMON_Q_BIN_EDGES[bin_index + 1]),
                    "q_units": "nm^-1",
                    "selected_as_sensitivity_peak": bin_index == common_peak_index,
                    "S_psi_psi_over_2_final_mean": _float(np.mean(bin_values)),
                    "S_psi_psi_over_2_final_sem": _float(_sem(bin_values)),
                    "S_psi_psi_over_2_early_mean": _float(np.mean(bin_early)),
                    "S_psi_psi_over_2_early_sem": _float(_sem(bin_early)),
                    "early_to_final_change_mean": _float(np.mean(bin_change)),
                    "early_to_final_change_sem": _float(_sem(bin_change)),
                    "final_window_slope_mean_per_step": _float(np.mean(bin_slope)),
                    "final_window_slope_sem_per_step": _float(_sem(bin_slope)),
                }
            )

    analysis_dir = root / "analysis"
    per_seed_path = analysis_dir / "per_seed_summary.csv"
    condition_path = analysis_dir / "condition_summary.csv"
    shell_seed_path = analysis_dir / "shell_spectra_by_seed.csv"
    shell_condition_path = analysis_dir / "shell_spectra_condition.csv"
    common_bin_condition_path = analysis_dir / "common_q_bin_sensitivity.csv"
    summary_path = analysis_dir / "summary.json"
    audit_path = analysis_dir / "analysis_manifest.json"

    _write_csv(per_seed_path, list(per_seed_rows[0]), per_seed_rows)
    _write_csv(condition_path, list(condition_rows[0]), condition_rows)
    _write_csv(shell_seed_path, list(shell_seed_rows[0]), shell_seed_rows)
    _write_csv(shell_condition_path, list(shell_condition_rows[0]), shell_condition_rows)
    _write_csv(
        common_bin_condition_path,
        list(common_bin_condition_rows[0]),
        common_bin_condition_rows,
    )
    summary = {
        "schema_version": 1,
        "campaign_id": manifest["campaign_id"],
        "campaign_design_sha256": manifest["design_sha256"],
        "n_planned_runs": len(plan),
        "n_complete_runs": len(run_records),
        "missing_runs": missing_runs,
        "primary_peak_metric": (
            "maximum of the across-seed mean final-five-frame exact-shell "
            "S_psi_psi^(N)(q)/2 spectrum"
        ),
        "per_seed_maxima_role": "secondary",
        "common_q_bin_edges_inverse_nm": COMMON_Q_BIN_EDGES.tolist(),
        "common_q_bin_candidate_count": len(COMMON_Q_BIN_EDGES) - 1,
        "diagnostic_metric": "S_AA(q)-S_AB(q), reported separately and not used for primary k*",
        "wavevector_units": "nm^-1",
        "uncertainty": "standard error across independent seeds (ddof=1)",
        "conditions": json_conditions,
    }
    atomic_json(summary_path, summary)

    output_paths = [
        per_seed_path,
        condition_path,
        shell_seed_path,
        shell_condition_path,
        common_bin_condition_path,
        summary_path,
    ]
    audit = {
        "schema_version": 1,
        "campaign_manifest": str(manifest_path.relative_to(root)),
        "campaign_manifest_sha256": sha256_file(manifest_path),
        "analysis_provenance": analysis_provenance,
        "analysis_provenance_sha256": sha256_bytes(canonical_json(analysis_provenance)),
        "source_files": source_files,
        "formulas": {
            "frame_reduction": "final five frames only define the estimator",
            "primary_peak": (
                "argmax exact shell of the across-seed mean final-five-frame "
                "S_psi_psi^(N)(q)/2 spectrum"
            ),
            "per_seed_peak": "secondary only",
            "condition_uncertainty": "seed SEM evaluated at the condition-selected shell",
            "common_q_bin_sensitivity": (
                "nine fixed physical-q bins with edges 0.15,0.30,...,1.50 nm^-1"
            ),
            "convergence": (
                "immediately preceding five-frame mean versus final five-frame mean; "
                "OLS slope across final five frames"
            ),
            "seed_sem": "sample standard deviation across seeds (ddof=1) divided by sqrt(n)",
            "diagnostic_only": "S_AA(q)-S_AB(q); never used to define the primary peak",
        },
        "outputs": {
            str(path.relative_to(root)): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in output_paths
        },
    }
    atomic_json(audit_path, audit)
    return {
        "per_seed_summary": per_seed_path,
        "condition_summary": condition_path,
        "shell_spectra_by_seed": shell_seed_path,
        "shell_spectra_condition": shell_condition_path,
        "common_q_bin_sensitivity": common_bin_condition_path,
        "summary_json": summary_path,
        "analysis_manifest": audit_path,
    }
