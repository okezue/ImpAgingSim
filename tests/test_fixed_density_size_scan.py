from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

import numpy as np

from melt.direct_structure import (
    DirectStructureFactorRecorder,
    direct_partial_structure_factors,
    early_and_final_snapshot_steps,
    final_snapshot_steps,
    reciprocal_modes,
)
from melt.fixed_density_analysis import COMMON_Q_BIN_EDGES, aggregate_campaign
from melt.fixed_density_size_scan import (
    DEFAULT_SEEDS,
    SCHEMA_VERSION,
    build_plan,
    campaign_design,
    canonical_json,
    completion_valid,
    ensure_campaign,
    expected_output_names,
    main,
    parse_args,
    run_spec,
    sha256_file,
    sha256_bytes,
    simulation_settings,
    validate_args,
)
from melt.sequences import generate_per_chain_exact_total
from melt.run import make_rng_streams


TEST_COMMIT = "a" * 40


def _clean_provenance() -> dict:
    return {
        "git_commit": TEST_COMMIT,
        "git_dirty": False,
        "git_diff_sha256": None,
    }


class TestFixedDensityPlan(unittest.TestCase):
    def test_default_plan_has_30_unique_fixed_density_runs(self):
        plan = build_plan()
        self.assertEqual(len(plan), 30)
        self.assertEqual(len({spec.run_id for spec in plan}), 30)
        self.assertEqual({spec.n_chains for spec in plan}, {144, 288, 576})
        self.assertEqual({spec.kappa for spec in plan}, {0.0, 1.0})
        self.assertEqual({spec.seed for spec in plan}, set(DEFAULT_SEEDS))
        density = 144 * 40 / 22.0**3
        for spec in plan:
            self.assertAlmostEqual(spec.box_size, 22.0 * (spec.n_chains / 144) ** (1 / 3), places=13)
            self.assertAlmostEqual(spec.bead_density, density, places=13)

    def test_grid_scaling_keeps_mesh_spacing_nearly_fixed(self):
        plan = build_plan(kappas=[0.0], seeds=[1])
        spacings = np.asarray([spec.grid_spacing for spec in plan])
        self.assertLess(float(np.ptp(spacings)), 0.003)
        self.assertLess(float(np.max(np.abs(spacings / spacings[0] - 1.0))), 0.007)
        self.assertEqual([spec.grid_size for spec in plan], [56, 71, 89])

    def test_common_physical_q_bins_are_populated_at_every_size(self):
        for spec in build_plan(kappas=[0.0], seeds=[1]):
            q = reciprocal_modes(spec.box_size, 1.5).q
            counts = []
            for index, (lower, upper) in enumerate(
                zip(COMMON_Q_BIN_EDGES[:-1], COMMON_Q_BIN_EDGES[1:])
            ):
                mask = (q >= lower) & (
                    (q <= upper) if index == len(COMMON_Q_BIN_EDGES) - 2 else (q < upper)
                )
                counts.append(int(np.sum(mask)))
            self.assertEqual(len(counts), 9)
            self.assertTrue(all(count > 0 for count in counts))

    def test_fixed_size_cli_rejects_asymmetric_composition_or_like_interactions(self):
        with self.assertRaisesRegex(ValueError, "requires f_A=0.5"):
            validate_args(parse_args(["--dry-run", "--f-A", "0.4"]))
        with self.assertRaisesRegex(ValueError, "requires lj_eps_AA=lj_eps_BB"):
            validate_args(parse_args(["--dry-run", "--lj-eps-BB", "0.9"]))

    def test_fixed_size_cli_rejects_partial_snapshot_interval(self):
        with self.assertRaisesRegex(ValueError, "exactly divisible"):
            validate_args(
                parse_args(["--dry-run", "--n-steps", "250001", "--snapshot-interval", "2000"])
            )

    def test_dry_run_enumerates_without_writing_or_openmm(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "campaigns"
            stream = io.StringIO()
            with redirect_stdout(stream):
                status = main(["--dry-run", "--out", str(out)])
            payload = json.loads(stream.getvalue())
            self.assertEqual(status, 0)
            self.assertEqual(payload["n_runs"], 30)
            self.assertFalse(out.exists())

    def test_manifest_only_writes_plan_without_running(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "campaigns"
            with patch(
                "melt.fixed_density_size_scan.collect_provenance",
                return_value=_clean_provenance(),
            ), redirect_stdout(io.StringIO()):
                status = main(["--manifest-only", "--out", str(out)])
            manifest = out / "fixed_density_pi099" / "manifest.json"
            self.assertEqual(status, 0)
            self.assertTrue(manifest.is_file())
            data = json.loads(manifest.read_text())
            self.assertEqual(data["design"]["n_runs"], 30)
            self.assertEqual(len(data["runs"]), 30)
            self.assertFalse(any((manifest.parent / "runs").iterdir()))

    def test_mocked_run_is_atomically_completed_then_checksum_cached(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = parse_args(
                [
                    "--run",
                    "--out",
                    temporary,
                    "--campaign-id",
                    "test",
                    "--sizes",
                    "144",
                    "--kappas",
                    "0",
                    "--seeds",
                    "7",
                    "--n-steps",
                    "1000",
                    "--snapshot-interval",
                    "100",
                ]
            )
            validate_args(args)
            plan = build_plan(
                sizes=args.sizes,
                kappas=args.kappas,
                seeds=args.seeds,
                chain_length=args.chain_length,
                reference_chains=args.reference_chains,
                reference_box_size=args.reference_box_size,
                reference_grid_size=args.reference_grid_size,
                pi=args.pi,
            )
            repo_root = Path(__file__).resolve().parents[1]
            campaign = Path(temporary) / "test"
            with patch(
                "melt.fixed_density_size_scan.collect_provenance",
                return_value=_clean_provenance(),
            ):
                ensure_campaign(campaign, campaign_design(args, plan), repo_root)

            def fake_execute(config):
                run_dir = Path(config.out) / config.run_id
                (run_dir / "meta.json").write_text("{}\n")
                (run_dir / "snapshots.csv").write_text("step\n1000\n")
                np.savez_compressed(run_dir / "structure_factor.npz", x=np.asarray([1]))
                np.savez_compressed(run_dir / "direct_structure_factor.npz", x=np.asarray([1]))
                return str(run_dir)

            with patch("melt.run.execute", side_effect=fake_execute) as execute_mock, patch(
                "melt.fixed_density_size_scan.collect_provenance",
                return_value=_clean_provenance(),
            ):
                status1, output1 = run_spec(plan[0], args, campaign, repo_root)
                status2, output2 = run_spec(plan[0], args, campaign, repo_root)
            self.assertEqual(status1, "completed")
            self.assertEqual(status2, "cached")
            self.assertEqual(output1, output2)
            self.assertTrue(completion_valid(output1))
            self.assertEqual(execute_mock.call_count, 1)

    def test_manifest_creation_rejects_dirty_tree(self):
        dirty = {**_clean_provenance(), "git_dirty": True, "git_diff_sha256": "b" * 64}
        with tempfile.TemporaryDirectory() as temporary, patch(
            "melt.fixed_density_size_scan.collect_provenance", return_value=dirty
        ):
            with self.assertRaisesRegex(RuntimeError, "git tree is dirty"):
                main(["--manifest-only", "--out", temporary])

    def test_existing_campaign_rejects_commit_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = parse_args(["--manifest-only", "--out", temporary, "--campaign-id", "pin"])
            plan = build_plan()
            campaign = Path(temporary) / "pin"
            repo_root = Path(__file__).resolve().parents[1]
            with patch(
                "melt.fixed_density_size_scan.collect_provenance",
                return_value=_clean_provenance(),
            ):
                ensure_campaign(campaign, campaign_design(args, plan), repo_root)
            other = {**_clean_provenance(), "git_commit": "b" * 40}
            with patch(
                "melt.fixed_density_size_scan.collect_provenance", return_value=other
            ), self.assertRaisesRegex(RuntimeError, "differs from campaign commit"):
                ensure_campaign(campaign, campaign_design(args, plan), repo_root)


class TestDirectStructureFactor(unittest.TestCase):
    def test_reciprocal_first_shell_has_six_modes(self):
        modes = reciprocal_modes(22.0, 1.5)
        self.assertEqual(int(modes.shell_n2[0]), 1)
        self.assertEqual(int(modes.shell_degeneracy[0]), 6)
        self.assertAlmostEqual(float(modes.shell_q[0]), 2.0 * np.pi / 22.0)
        self.assertTrue(np.all(modes.q <= 1.5 + 1e-14))

    def test_direct_partial_normalization_matches_two_bead_formula(self):
        positions = np.asarray([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
        types = np.asarray([1, 0], dtype=np.int8)
        q = np.asarray([[np.pi, 0.0, 0.0], [2.0 * np.pi, 0.0, 0.0]])
        spectra = direct_partial_structure_factors(positions, types, q, chunk_size=1)
        np.testing.assert_allclose(spectra["S_AA"], [0.5, 0.5], atol=1e-15)
        np.testing.assert_allclose(spectra["S_BB"], [0.5, 0.5], atol=1e-15)
        np.testing.assert_allclose(spectra["S_AB"], [0.0, -0.5], atol=1e-15)
        np.testing.assert_allclose(spectra["S_CC"], [1.0, 2.0], atol=1e-15)
        np.testing.assert_allclose(spectra["S_psi_psi_over_2"], [0.5, 1.0], atol=1e-15)

    def test_primary_channel_is_label_swap_invariant_but_A_diagnostic_is_not(self):
        positions = np.asarray([[0.0, 0.0, 0.0], [0.23, 0.0, 0.0], [0.81, 0.0, 0.0]])
        types = np.asarray([1, 1, 0], dtype=np.int8)
        q = np.asarray([[2.7, 0.0, 0.0], [4.1, 0.0, 0.0]])
        original = direct_partial_structure_factors(positions, types, q)
        swapped = direct_partial_structure_factors(positions, 1 - types, q)
        np.testing.assert_allclose(
            original["S_psi_psi_over_2"], swapped["S_psi_psi_over_2"], rtol=1e-14
        )
        original_diagnostic = original["S_AA"] - original["S_AB"]
        swapped_diagnostic = swapped["S_AA"] - swapped["S_AB"]
        self.assertFalse(np.allclose(original_diagnostic, swapped_diagnostic, rtol=1e-12))

    def test_recorder_selects_only_requested_final_steps_and_saves_shells(self):
        rng = np.random.default_rng(5)
        positions = rng.uniform(0.0, 10.0, size=(20, 3))
        types = np.asarray([0, 1] * 10, dtype=np.int8)
        steps = final_snapshot_steps(1000, 100, 2)
        recorder = DirectStructureFactorRecorder(10.0, 1.5, steps, chunk_size=4)
        self.assertFalse(recorder.observe(800, positions, types))
        self.assertTrue(recorder.observe(900, positions, types))
        self.assertTrue(recorder.observe(1000, positions, types))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "direct.npz"
            recorder.save(str(output))
            with np.load(output, allow_pickle=False) as data:
                np.testing.assert_array_equal(data["steps"], [900, 1000])
                self.assertEqual(data["S_AA"].shape[0], 2)
                self.assertEqual(data["S_AA_shell"].shape, (2, len(data["shell_q"])))

    def test_recorder_rejects_unsorted_target_steps(self):
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            DirectStructureFactorRecorder(
                10.0,
                1.5,
                np.asarray([1000, 900]),
                window_labels=np.asarray(["final", "early"]),
            )

    def test_final_snapshot_steps_rejects_more_frames_than_exist(self):
        with self.assertRaisesRegex(ValueError, "records only 2 snapshots"):
            final_snapshot_steps(200, 100, 3)

    def test_production_schedule_selects_exact_final_five_frames(self):
        np.testing.assert_array_equal(
            final_snapshot_steps(250000, 2000, 5),
            [242000, 244000, 246000, 248000, 250000],
        )
        early, final = early_and_final_snapshot_steps(250000, 2000, 5, 5)
        np.testing.assert_array_equal(early, [232000, 234000, 236000, 238000, 240000])
        np.testing.assert_array_equal(final, [242000, 244000, 246000, 248000, 250000])


class TestExactGlobalComposition(unittest.TestCase):
    def test_rejection_sampler_is_exact_and_deterministic(self):
        kwargs = ("correlated", 12, 10, 0.5, 4, 0.0, 0.99)
        first, first_attempts = generate_per_chain_exact_total(
            *kwargs, np.random.default_rng(123), 10000
        )
        second, second_attempts = generate_per_chain_exact_total(
            *kwargs, np.random.default_rng(123), 10000
        )
        np.testing.assert_array_equal(first, second)
        self.assertEqual(first_attempts, second_attempts)
        self.assertEqual(float(np.mean(first)), 0.5)

    def test_retained_chains_are_not_mechanically_complement_paired(self):
        sequence, _ = generate_per_chain_exact_total(
            "correlated", 20, 20, 0.5, 4, 0.0, 0.99, np.random.default_rng(123), 10000
        )
        chains = sequence.reshape(20, 20)
        complement_matches = [
            any(np.array_equal(1 - chain, other) for other in chains) for chain in chains
        ]
        self.assertFalse(any(complement_matches))

    def test_impossible_target_and_attempt_exhaustion_raise(self):
        with self.assertRaisesRegex(ValueError, "integer target"):
            generate_per_chain_exact_total(
                "correlated", 3, 3, 0.5, 4, 0.0, 0.99, np.random.default_rng(1), 10
            )
        with self.assertRaisesRegex(RuntimeError, "after 1 independent attempts"):
            generate_per_chain_exact_total(
                "correlated", 12, 10, 0.5, 4, 0.0, 0.99, np.random.default_rng(123), 1
            )

    def test_split_placement_stream_is_independent_of_sequence_attempt_consumption(self):
        sequence_rng_1, placement_rng_1, metadata_1 = make_rng_streams(17, split=True)
        sequence_rng_1.random(10)
        placement_draw_1 = placement_rng_1.random(12)
        sequence_rng_2, placement_rng_2, metadata_2 = make_rng_streams(17, split=True)
        sequence_rng_2.random(10000)
        placement_draw_2 = placement_rng_2.random(12)
        np.testing.assert_array_equal(placement_draw_1, placement_draw_2)
        self.assertEqual(metadata_1, metadata_2)


def _synthetic_direct(
    path: Path,
    spec: dict,
    settings: dict,
    peak_shell: int,
    amplitude: float,
) -> None:
    early_steps, final_steps = early_and_final_snapshot_steps(
        settings["n_steps"],
        settings["snapshot_interval"],
        settings["direct_early_frames"],
        settings["direct_final_frames"],
    )
    steps = np.concatenate((early_steps, final_steps))
    windows = np.asarray(["early"] * len(early_steps) + ["final"] * len(final_steps))
    modes = reciprocal_modes(spec["box_size"], settings["direct_q_max"])
    shell_AA = np.ones(len(modes.shell_q), dtype=float)
    shell_AA[peak_shell] += amplitude
    mode_AA = shell_AA[modes.shell_index]
    S_AA = np.tile(mode_AA, (len(steps), 1))
    S_BB = S_AA.copy()
    S_AB = np.zeros_like(S_AA)
    S_CC = S_AA + S_BB - 2.0 * S_AB
    S_primary = 0.5 * S_CC
    S_AA_shell = np.tile(shell_AA, (len(steps), 1))
    S_BB_shell = S_AA_shell.copy()
    S_AB_shell = np.zeros_like(S_AA_shell)
    S_CC_shell = S_AA_shell + S_BB_shell - 2.0 * S_AB_shell
    S_primary_shell = 0.5 * S_CC_shell
    np.savez_compressed(
        path,
        steps=steps,
        window=windows,
        h=modes.h,
        q_vectors=modes.q_vectors,
        q=modes.q,
        n2=modes.n2,
        shell_n2=modes.shell_n2,
        shell_q=modes.shell_q,
        shell_index=modes.shell_index,
        shell_degeneracy=modes.shell_degeneracy,
        box_size_nm=np.asarray(spec["box_size"]),
        q_max_inverse_nm=np.asarray(settings["direct_q_max"]),
        S_AA=S_AA,
        S_BB=S_BB,
        S_AB=S_AB,
        S_CC=S_CC,
        S_psi_psi_over_2=S_primary,
        S_AA_shell=S_AA_shell,
        S_BB_shell=S_BB_shell,
        S_AB_shell=S_AB_shell,
        S_CC_shell=S_CC_shell,
        S_psi_psi_over_2_shell=S_primary_shell,
        normalization=np.asarray("S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total"),
        composition_channel=np.asarray("S_psi_psi^(N)=S_CC=S_AA+S_BB-2*S_AB"),
        primary_channel=np.asarray("S_psi_psi^(N)/2=S_CC/2"),
        q_units=np.asarray("nm^-1"),
    )


class TestFixedDensityAggregation(unittest.TestCase):
    def _campaign(self, root: Path) -> Path:
        campaign = root / "campaign"
        settings_args = parse_args(
            ["--dry-run", "--n-steps", "20000", "--snapshot-interval", "2000"]
        )
        settings = simulation_settings(settings_args)
        runs = []
        index = 0
        for seed, peak_shell, amplitude in (
            (11, 1, 2.0),
            (13, 1, 4.0),
        ):
            run_id = f"run_{seed}"
            spec = {
                "index": index,
                "run_id": run_id,
                "n_chains": 144,
                "chain_length": 40,
                "n_beads": 5760,
                "box_size": 22.0,
                "bead_density": 5760 / 22.0**3,
                "grid_size": 56,
                "grid_spacing": 22.0 / 56,
                "kappa": 1.0,
                "pi": 0.99,
                "seed": seed,
            }
            runs.append(spec)
            run_dir = campaign / "runs" / run_id
            run_dir.mkdir(parents=True)
            direct = run_dir / "direct_structure_factor.npz"
            _synthetic_direct(direct, spec, settings, peak_shell, amplitude)
            n_A = spec["n_beads"] // 2
            sequence = [1] * n_A + [0] * (spec["n_beads"] - n_A)
            meta = {
                "melt_params": {
                    "n_chains": spec["n_chains"],
                    "chain_length": spec["chain_length"],
                    "box_size": spec["box_size"],
                    "bond_k": settings["bond_k"],
                    "bond_r0": settings["bond_r0"],
                    "lj_eps_AA": settings["lj_eps_AA"],
                    "lj_eps_BB": settings["lj_eps_BB"],
                    "lj_eps_AB": settings["lj_eps_AB"],
                    "lj_eps_core": settings["lj_eps_core"],
                    "lj_sigma": settings["lj_sigma"],
                    "lj_cutoff": settings["lj_cutoff"],
                    "temperature": settings["T_quench"] * settings["lj_eps_AA"] / 0.00831446261815324,
                    "friction": settings["friction"],
                    "dt": settings["dt"],
                },
                "run_params": {
                    "n_steps": settings["n_steps"],
                    "equilibration_steps": settings["equilibration"],
                    "snapshot_interval": settings["snapshot_interval"],
                    "seed": seed,
                },
                "sequence_type": settings["sequence"],
                "sequence": sequence,
                "n_particles": spec["n_beads"],
                "T_equilibrate_star": settings["T_equilibrate"],
                "T_quench_star": settings["T_quench"],
                "T_equilibrate_kelvin": (
                    settings["T_equilibrate"] * settings["lj_eps_AA"] / 0.00831446261815324
                ),
                "T_quench_kelvin": (
                    settings["T_quench"] * settings["lj_eps_AA"] / 0.00831446261815324
                ),
                "eps_ref_kjmol": settings["lj_eps_AA"],
                "k_B_kjmolK": 0.00831446261815324,
                "random_streams": {
                    "method": "numpy_SeedSequence_spawn",
                    "root_seed": seed,
                    "sequence_spawn_key": [0],
                    "placement_spawn_key": [1],
                    "openmm_integrator_seed": seed,
                },
                "sequence_parameters": {
                    "f_A": settings["f_A"],
                    "kappa": spec["kappa"],
                    "pi": spec["pi"],
                    "block_length": settings["block_length"],
                },
                "density_analysis": {"enabled": True, "grid_size": spec["grid_size"]},
                "requested_openmm_platform": settings["platform"],
                "sequence_ensemble": {
                    "method": (
                        "independent_per_chain_rejection_conditioned_on_exact_global_count"
                    ),
                    "n_chains_drawn_independently_per_attempt": spec["n_chains"],
                    "acceptance_attempts": 3 if seed == 11 else 5,
                    "max_attempts": settings["sequence_balance_max_attempts"],
                    "target_A_count": n_A,
                    "realized_A_count": n_A,
                    "exact_global_f_A": 0.5,
                    "constructed_complement_pairs": False,
                },
            }
            (run_dir / "meta.json").write_text(json.dumps(meta))
            (run_dir / "planned_run.json").write_text(
                json.dumps({"run_spec": spec, "simulation_settings": settings})
            )
            (run_dir / "snapshots.csv").write_text("step\n20000\n")
            np.savez_compressed(run_dir / "structure_factor.npz", x=np.asarray([1]))
            required_names = expected_output_names(False)
            files = {
                name: {
                    "bytes": (run_dir / name).stat().st_size,
                    "sha256": sha256_file(run_dir / name),
                }
                for name in required_names
            }
            completion = {
                "schema_version": SCHEMA_VERSION,
                "run_spec": spec,
                "save_trajectory": False,
                "required_outputs": list(required_names),
                "provenance": _clean_provenance(),
                "files": files,
            }
            (run_dir / "completion.json").write_text(json.dumps(completion))
            index += 1
        design_payload = {
            "schema_version": SCHEMA_VERSION,
            "study": "synthetic fixed-density finite-size comparison",
            "campaign_id": "synthetic",
            "design": {"synthetic": True},
            "simulation_settings": settings,
            "direct_structure_factor": {
                "synthetic": True,
                "common_q_bin_edges_inverse_nm": COMMON_Q_BIN_EDGES.tolist(),
            },
            "runs": runs,
        }
        creation_provenance = _clean_provenance()
        manifest = {
            **design_payload,
            "design_sha256": sha256_bytes(canonical_json(design_payload)),
            "creation_provenance": creation_provenance,
            "creation_provenance_sha256": sha256_bytes(canonical_json(creation_provenance)),
        }
        campaign.mkdir(exist_ok=True)
        (campaign / "manifest.json").write_text(json.dumps(manifest))
        return campaign

    def _first_run(self, campaign: Path) -> tuple[Path, dict]:
        manifest = json.loads((campaign / "manifest.json").read_text())
        spec = manifest["runs"][0]
        return campaign / "runs" / spec["run_id"], spec

    def _refresh_completion_file(self, run_dir: Path, name: str) -> None:
        completion_path = run_dir / "completion.json"
        completion = json.loads(completion_path.read_text())
        path = run_dir / name
        completion["files"][name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        completion_path.write_text(json.dumps(completion))

    def test_completion_rejects_empty_and_missing_required_file_maps(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            run_dir, spec = self._first_run(campaign)
            completion_path = run_dir / "completion.json"
            original = json.loads(completion_path.read_text())
            self.assertTrue(completion_valid(run_dir, spec, False, TEST_COMMIT))
            empty = {**original, "files": {}}
            completion_path.write_text(json.dumps(empty))
            self.assertFalse(completion_valid(run_dir, spec, False, TEST_COMMIT))
            missing = json.loads(json.dumps(original))
            missing["files"].pop("structure_factor.npz")
            completion_path.write_text(json.dumps(missing))
            self.assertFalse(completion_valid(run_dir, spec, False, TEST_COMMIT))
            trajectory_claim = json.loads(json.dumps(original))
            trajectory_claim["save_trajectory"] = True
            trajectory_claim["required_outputs"] = list(expected_output_names(True))
            completion_path.write_text(json.dumps(trajectory_claim))
            self.assertFalse(completion_valid(run_dir, spec, True, TEST_COMMIT))

    def test_completion_rejects_schema_and_checksum_corruption(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            run_dir, spec = self._first_run(campaign)
            completion_path = run_dir / "completion.json"
            original = json.loads(completion_path.read_text())
            wrong_schema = {**original, "schema_version": SCHEMA_VERSION + 1}
            completion_path.write_text(json.dumps(wrong_schema))
            self.assertFalse(completion_valid(run_dir, spec, False, TEST_COMMIT))
            completion_path.write_text(json.dumps(original))
            with (run_dir / "meta.json").open("a") as handle:
                handle.write("corrupt")
            self.assertFalse(completion_valid(run_dir, spec, False, TEST_COMMIT))

    def test_aggregation_writes_seed_condition_shell_and_json_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ):
                outputs = aggregate_campaign(campaign)
            for path in outputs.values():
                self.assertTrue(path.is_file())
            with outputs["per_seed_summary"].open() as handle:
                seed_rows = list(csv.DictReader(handle))
            self.assertEqual(len(seed_rows), 2)
            self.assertAlmostEqual(
                float(seed_rows[0]["secondary_per_seed_max_k_star"]),
                np.sqrt(2) * 2.0 * np.pi / 22.0,
            )
            self.assertEqual(float(seed_rows[0]["realized_f_A"]), 0.5)
            with outputs["condition_summary"].open() as handle:
                condition = next(csv.DictReader(handle))
            self.assertAlmostEqual(
                float(condition["primary_amplitude_seed_mean_at_selected_shell"]), 4.0
            )
            self.assertAlmostEqual(
                float(condition["primary_amplitude_seed_sem_at_selected_shell"]), 1.0
            )
            self.assertEqual(float(condition["realized_f_A_mean"]), 0.5)
            self.assertEqual(float(condition["realized_f_A_sem"]), 0.0)
            self.assertEqual(int(condition["common_bin_candidate_count"]), 9)
            self.assertEqual(
                float(condition["primary_early_to_final_change_seed_mean_at_selected_shell"]),
                0.0,
            )
            self.assertEqual(
                float(condition["primary_final_window_slope_seed_mean_per_step_at_selected_shell"]),
                0.0,
            )
            summary = json.loads(outputs["summary_json"].read_text())
            self.assertIn("across-seed mean", summary["primary_peak_metric"])
            self.assertEqual(summary["n_complete_runs"], 2)
            self.assertEqual(summary["missing_runs"], [])

    def test_aggregation_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ):
                first = aggregate_campaign(campaign)
            hashes_before = {name: sha256_file(path) for name, path in first.items()}
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ):
                second = aggregate_campaign(campaign)
            hashes_after = {name: sha256_file(path) for name, path in second.items()}
            self.assertEqual(hashes_before, hashes_after)

    def test_aggregation_rejects_manifest_design_tamper(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            manifest_path = campaign / "manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["simulation_settings"]["direct_q_max"] = 1.4
            manifest_path.write_text(json.dumps(manifest))
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ), self.assertRaisesRegex(ValueError, "failed its design hash"):
                aggregate_campaign(campaign)

    def test_aggregation_rejects_planned_settings_mismatch_after_rehash(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            run_dir, _ = self._first_run(campaign)
            planned_path = run_dir / "planned_run.json"
            planned = json.loads(planned_path.read_text())
            planned["simulation_settings"]["friction"] = 9.0
            planned_path.write_text(json.dumps(planned))
            self._refresh_completion_file(run_dir, planned_path.name)
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ), self.assertRaisesRegex(ValueError, "simulation_settings do not match"):
                aggregate_campaign(campaign)

    def test_aggregation_rejects_meta_force_setting_mismatch_after_rehash(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            run_dir, _ = self._first_run(campaign)
            meta_path = run_dir / "meta.json"
            meta = json.loads(meta_path.read_text())
            meta["melt_params"]["bond_k"] = 999.0
            meta_path.write_text(json.dumps(meta))
            self._refresh_completion_file(run_dir, meta_path.name)
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ), self.assertRaisesRegex(ValueError, "bond_k does not match"):
                aggregate_campaign(campaign)

    def test_aggregation_rejects_current_or_run_commit_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            other = {**_clean_provenance(), "git_commit": "b" * 40}
            with patch(
                "melt.fixed_density_analysis.collect_provenance", return_value=other
            ), self.assertRaisesRegex(RuntimeError, "differs from campaign commit"):
                aggregate_campaign(campaign)

        with tempfile.TemporaryDirectory() as temporary:
            campaign = self._campaign(Path(temporary))
            run_dir, _ = self._first_run(campaign)
            completion_path = run_dir / "completion.json"
            completion = json.loads(completion_path.read_text())
            completion["provenance"] = {**_clean_provenance(), "git_commit": "b" * 40}
            completion_path.write_text(json.dumps(completion))
            with patch(
                "melt.fixed_density_analysis.collect_provenance",
                return_value=_clean_provenance(),
            ), self.assertRaisesRegex(RuntimeError, "planned runs are missing or invalid"):
                aggregate_campaign(campaign)

    def test_aggregation_rejects_step_and_reciprocal_metadata_tamper(self):
        for field in ("steps", "q_vectors"):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                campaign = self._campaign(Path(temporary))
                run_dir, _ = self._first_run(campaign)
                direct_path = run_dir / "direct_structure_factor.npz"
                with np.load(direct_path, allow_pickle=False) as archive:
                    payload = {name: np.asarray(archive[name]) for name in archive.files}
                payload[field] = payload[field].copy()
                payload[field].flat[0] += 1
                np.savez_compressed(direct_path, **payload)
                self._refresh_completion_file(run_dir, direct_path.name)
                with patch(
                    "melt.fixed_density_analysis.collect_provenance",
                    return_value=_clean_provenance(),
                ), self.assertRaises(AssertionError):
                    aggregate_campaign(campaign)


if __name__ == "__main__":
    unittest.main()
