from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from melt.campaign import completion_valid, design_hash, select_shard
from melt.epsab_scan import (
    DEFAULT_EPS_ABS,
    DESIGN_KEYS,
    REQUIRED_OUTPUTS,
    build_plan,
    campaign_design,
    make_run_id,
    main,
    parse_args,
)

try:
    import openmm  # noqa: F401

    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
needs_openmm = pytest.mark.skipif(not HAS_OPENMM, reason="openmm not installed")


class TestPlan:
    def test_default_grid_spans_chi_zero_to_production(self):
        assert DEFAULT_EPS_ABS[0] == 1.0
        assert abs(DEFAULT_EPS_ABS[-1] - 0.1) < 1e-9
        assert len(DEFAULT_EPS_ABS) == 37
        assert np.allclose(np.diff(DEFAULT_EPS_ABS), -0.025)

    def test_plan_size_and_box_scaling(self):
        plan = build_plan(eps_ABs=[1.0, 0.5], seeds=[1, 2, 3], sizes=[144, 288])
        assert len(plan) == 12
        assert [s.index for s in plan] == list(range(12))
        by_size = {s.n_chains: s for s in plan}
        assert abs(by_size[144].box_size - 22.0) < 1e-12
        assert abs(by_size[288].box_size - 22.0 * 2 ** (1 / 3)) < 1e-9
        assert abs(by_size[144].bead_density - by_size[288].bead_density) < 1e-12
        assert len({s.run_id for s in plan}) == 12

    def test_run_id_format(self):
        assert make_run_id(144, 40, 0.5, 0.99, 0.975, 3) == "M0144_N40_kappa0p5_pi0p990_epsAB0p975_seed3"
        assert make_run_id(144, 40, 0.5, 0.99, 1.0, 1) == "M0144_N40_kappa0p5_pi0p990_epsAB1_seed1"

    def test_plan_rejects_bad_inputs(self):
        with pytest.raises(ValueError):
            build_plan(eps_ABs=[0.5, 0.5])
        with pytest.raises(ValueError):
            build_plan(sizes=[143])
        with pytest.raises(ValueError):
            build_plan(kappa=1.5)

    def test_design_hash_changes_with_grid(self):
        a = parse_args(["--dry-run", "--eps-ABs", "1.0", "0.5"])
        b = parse_args(["--dry-run", "--eps-ABs", "1.0", "0.4"])
        ha = design_hash(campaign_design(a, build_plan(a.eps_ABs, a.seeds, a.sizes)), DESIGN_KEYS)
        hb = design_hash(campaign_design(b, build_plan(b.eps_ABs, b.seeds, b.sizes)), DESIGN_KEYS)
        assert ha != hb
        # The platform is an execution detail and must not enter the design hash.
        c = parse_args(["--dry-run", "--eps-ABs", "1.0", "0.5", "--platform", "CUDA"])
        hc = design_hash(campaign_design(c, build_plan(c.eps_ABs, c.seeds, c.sizes)), DESIGN_KEYS)
        assert hc == ha

    def test_settings_pin_the_production_protocol(self):
        a = parse_args(["--dry-run"])
        settings = campaign_design(a, build_plan())["simulation_settings"]
        assert settings["n_steps"] == 250000 and settings["snapshot_interval"] == 2000
        assert settings["equilibration"] == 30000 and settings["T_quench"] == 0.7
        assert settings["lj_eps_AA"] == 1.0 and settings["lj_eps_BB"] == 1.0
        assert settings["record_modes"] is True and settings["mode_interval"] == 200
        assert settings["exact_global_composition"] is True


class TestParallelLauncher:
    def test_cuda_devices_round_robin(self, monkeypatch):
        from melt import epsab_scan

        seen = []

        class Proc:
            returncode = 0
            stdout = "ok"

        def fake_run(cmd, env, **kwargs):
            seen.append((int(cmd[cmd.index("--run-index") + 1]), env.get("CUDA_VISIBLE_DEVICES")))
            return Proc()

        monkeypatch.setattr(epsab_scan.subprocess, "run", fake_run)
        args = parse_args(["--run", "--parallel", "2", "--cuda-devices", "0,1,2"])
        assert epsab_scan.run_parallel(args, [10, 11, 12, 13]) == 0
        assert sorted(seen) == [(10, "0"), (11, "1"), (12, "2"), (13, "0")]

    def test_no_cuda_devices_leaves_env_alone(self, monkeypatch):
        from melt import epsab_scan

        envs = []

        class Proc:
            returncode = 0
            stdout = ""

        monkeypatch.setattr(epsab_scan.subprocess, "run", lambda cmd, env, **k: envs.append(env) or Proc())
        monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
        args = parse_args(["--run", "--parallel", "2", "--threads-per-run", "1"])
        assert epsab_scan.run_parallel(args, [0]) == 0
        assert "CUDA_VISIBLE_DEVICES" not in envs[0] and envs[0]["OPENMM_CPU_THREADS"] == "1"


class TestSharding:
    def test_shards_partition_indices(self):
        indices = list(range(10))
        shards = [select_shard(indices, k, 3) for k in range(3)]
        assert sorted(sum(shards, [])) == indices
        assert shards[0] == [0, 3, 6, 9]

    def test_shard_validation(self):
        with pytest.raises(ValueError):
            select_shard([0, 1], 0, None)
        with pytest.raises(ValueError):
            select_shard([0, 1], 3, 3)


SMOKE_ARGS = [
    "--eps-ABs", "1.0", "0.1", "--seeds", "1", "--sizes", "12", "--chain-length", "8",
    "--equilibration", "100", "--n-steps", "800", "--snapshot-interval", "400",
    "--mode-interval", "100", "--mode-q-max", "2.5", "--platform", "Reference", "--allow-dirty",
]


@needs_openmm
class TestCampaignEndToEnd:
    def test_manifest_run_cache_status_analyze(self, tmp_path):
        out = str(tmp_path / "camp")
        common = ["--out", out, "--campaign-id", "t", *SMOKE_ARGS]
        assert main(["--manifest-only", *common]) == 0
        manifest = json.loads((Path(out) / "t" / "manifest.json").read_text())
        assert manifest["design"]["n_runs"] == 2
        assert manifest["require_clean_git"] is False

        assert main(["--run", *common]) == 0
        runs = Path(out) / "t" / "runs"
        for spec in manifest["runs"]:
            run_dir = runs / spec["run_id"]
            assert completion_valid(run_dir, spec, REQUIRED_OUTPUTS, None)
            assert (run_dir / "mode_amplitudes.npz").is_file()

        # Corrupt one output: the run must be detected as incomplete and be redone on request.
        victim = runs / manifest["runs"][0]["run_id"]
        (victim / "snapshots.csv").write_text("corrupt\n")
        assert not completion_valid(victim, manifest["runs"][0], REQUIRED_OUTPUTS, None)
        with pytest.raises(RuntimeError):
            main(["--run", "--run-index", "0", *common])
        assert main(["--run", "--run-index", "0", "--recover-interrupted", *common]) == 0
        assert completion_valid(victim, manifest["runs"][0], REQUIRED_OUTPUTS, None)
        assert any((Path(out) / "t" / ".interrupted").iterdir())

        assert main(["--status", *common]) == 0
        assert main(["--analyze", *common]) == 0
        analysis = Path(out) / "t" / "analysis"
        assert (analysis / "per_run.csv").is_file()
        assert (analysis / "per_condition.csv").is_file()
        assert (analysis / "transition_summary.json").is_file()
        assert (analysis / "figures" / "S_peak_vs_delta_eps.png").is_file()
        summary = json.loads((analysis / "transition_summary.json").read_text())
        assert summary["n_runs_analyzed"] == 2
        assert "12" in summary["rpa"] and "12" in summary["transition"]

    def test_changed_design_needs_new_campaign_id(self, tmp_path):
        out = str(tmp_path / "camp")
        assert main(["--manifest-only", "--out", out, "--campaign-id", "t", *SMOKE_ARGS]) == 0
        changed = [a if a != "0.1" else "0.2" for a in SMOKE_ARGS]
        with pytest.raises(RuntimeError):
            main(["--manifest-only", "--out", out, "--campaign-id", "t", *changed])

    def test_shards_cover_the_plan(self, tmp_path):
        out = str(tmp_path / "camp")
        common = ["--out", out, "--campaign-id", "t", *SMOKE_ARGS]
        assert main(["--manifest-only", *common]) == 0
        assert main(["--run", "--shard-index", "0", "--shard-count", "2", *common]) == 0
        assert main(["--run", "--shard-index", "1", "--shard-count", "2", *common]) == 0
        status_manifest = json.loads((Path(out) / "t" / "manifest.json").read_text())
        for spec in status_manifest["runs"]:
            assert completion_valid(Path(out) / "t" / "runs" / spec["run_id"], spec, REQUIRED_OUTPUTS, None)
