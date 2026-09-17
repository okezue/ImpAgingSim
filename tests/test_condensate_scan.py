from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from melt.campaign import completion_valid, design_hash
from melt.condensate_scan import DESIGN_KEYS, build_plan, campaign_design, main, parse_args, required_outputs

try:
    import openmm  # noqa: F401

    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
needs_openmm = pytest.mark.skipif(not HAS_OPENMM, reason="openmm not installed")


class TestPlan:
    def test_protein_plan_and_box(self):
        args = parse_args(["--dry-run", "--eps-BBs", "1.0", "2.0", "--kappas", "0", "1", "--f-As", "0.7", "--seeds", "1", "2",
                           "--bead-density", "0.2"])
        plan = build_plan(args)
        assert len(plan) == 8
        assert abs(plan[0].box_size - (5760 / 0.2) ** (1 / 3)) < 1e-9
        assert all(not s.marks_dynamic and s.k_off == 0.0 and s.k_fb == 0.0 for s in plan)
        assert plan[0].run_id == "epsBB1_kappa0_fA0p7_seed1"
        assert len({s.run_id for s in plan}) == 8

    def test_chromatin_plan_auto_k_on(self):
        args = parse_args(["--dry-run", "--marks-dynamic", "--eps-BBs", "1.5", "--kappas", "0.5", "--f-As", "0.7",
                           "--k-offs", "0.01", "--k-fbs", "0", "0.1", "--seeds", "1"])
        plan = build_plan(args)
        assert len(plan) == 2
        # Without feedback the steady B fraction is k_on/(k_on+k_off) = 0.3.
        for s in plan:
            assert abs(s.k_on / (s.k_on + s.k_off) - 0.3) < 1e-9
        assert "koff0p01_kfb0p1" in plan[1].run_id
        assert required_outputs(True, False)[-1] == "marks.npz"

    def test_design_hash_depends_on_mode_and_grid(self):
        a = parse_args(["--dry-run", "--eps-BBs", "1.0", "--kappas", "0.5", "--seeds", "1"])
        b = parse_args(["--dry-run", "--marks-dynamic", "--eps-BBs", "1.0", "--kappas", "0.5", "--seeds", "1"])
        ha = design_hash(campaign_design(a, build_plan(a)), DESIGN_KEYS)
        hb = design_hash(campaign_design(b, build_plan(b)), DESIGN_KEYS)
        assert ha != hb
        c = parse_args(["--dry-run", "--eps-BBs", "1.0", "--kappas", "0.5", "--seeds", "1", "--platform", "CUDA"])
        assert design_hash(campaign_design(c, build_plan(c)), DESIGN_KEYS) == ha

    def test_validation(self):
        from melt.condensate_scan import validate_args

        with pytest.raises(ValueError):
            validate_args(parse_args(["--dry-run", "--f-As", "0.7", "--n-chains", "12", "--chain-length", "8"]))
        with pytest.raises(ValueError):
            validate_args(parse_args(["--dry-run", "--condensation-interval", "3000", "--snapshot-interval", "2000"]))
        with pytest.raises(ValueError):
            validate_args(parse_args(["--dry-run", "--bead-density", "1.5"]))

    def test_energetics_are_asymmetric(self):
        args = parse_args(["--dry-run"])
        settings = campaign_design(args, build_plan(args))["simulation_settings"]
        assert settings["lj_eps_AA"] == 0.0 and settings["lj_eps_AB"] == 0.0 and settings["lj_eps_core"] == 1.0
        assert settings["record_modes"] is True and settings["condensation_interval"] > 0


SMOKE = ["--eps-BBs", "1.0", "2.0", "--kappas", "1.0", "--f-As", "0.75", "--seeds", "1", "--n-chains", "12", "--chain-length", "8",
         "--bead-density", "0.2", "--equilibration", "100", "--n-steps", "800", "--snapshot-interval", "400",
         "--mode-interval", "100", "--condensation-interval", "400", "--mark-interval", "100", "--mode-q-max", "2.5",
         "--platform", "Reference", "--allow-dirty"]


@needs_openmm
class TestEndToEnd:
    def test_protein_campaign(self, tmp_path):
        out = str(tmp_path / "c")
        common = ["--out", out, "--campaign-id", "p", *SMOKE]
        assert main(["--manifest-only", *common]) == 0
        assert main(["--run", *common]) == 0
        manifest = json.loads((Path(out) / "p" / "manifest.json").read_text())
        for spec in manifest["runs"]:
            assert completion_valid(Path(out) / "p" / "runs" / spec["run_id"], spec, required_outputs(False, False), None)
        assert main(["--analyze", *common]) == 0
        rows = list(csv.DictReader(open(Path(out) / "p" / "analysis" / "per_condition.csv")))
        assert len(rows) == 2
        assert all(r["largest_B_cluster_fraction_mean"] != "" for r in rows)
        assert (Path(out) / "p" / "analysis" / "figures" / "heat_largest_cluster.png").is_file()

    def test_chromatin_campaign(self, tmp_path):
        out = str(tmp_path / "c")
        common = ["--out", out, "--campaign-id", "m", "--marks-dynamic", "--k-offs", "0.5", "--k-fbs", "0", "2.0", *SMOKE]
        assert main(["--manifest-only", *common]) == 0
        assert main(["--run", *common]) == 0
        manifest = json.loads((Path(out) / "m" / "manifest.json").read_text())
        assert len(manifest["runs"]) == 4
        for spec in manifest["runs"]:
            run_dir = Path(out) / "m" / "runs" / spec["run_id"]
            assert completion_valid(run_dir, spec, required_outputs(True, False), None)
            with np.load(run_dir / "marks.npz") as d:
                assert d["types"].shape[0] == 8
        assert main(["--analyze", *common]) == 0
        rows = list(csv.DictReader(open(Path(out) / "m" / "analysis" / "per_condition.csv")))
        assert len(rows) == 4
        assert all(r["mark_memory_time_mean"] != "" or r["f_B_mean"] != "" for r in rows)
        assert (Path(out) / "m" / "analysis" / "figures" / "heat_f_B.png").is_file()
