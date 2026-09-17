from __future__ import annotations

import csv
import json
import os

import numpy as np
import pytest

from melt.clusters import cluster_labels, condensation_summary, contact_pairs
from melt.marks import (
    TYPE_A,
    TYPE_B,
    MarkDynamics,
    MarkRecorder,
    hill,
    local_B_counts,
    mark_autocorrelation,
    mean_field_rates,
    persistence_fraction,
    step_marks,
)

try:
    import openmm  # noqa: F401

    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
needs_openmm = pytest.mark.skipif(not HAS_OPENMM, reason="openmm not installed")


class TestLocalCounts:
    def test_counts_exclude_self_and_wrap_periodically(self):
        L = 10.0
        pos = np.array([[0.2, 5.0, 5.0], [9.9, 5.0, 5.0], [5.0, 5.0, 5.0], [0.5, 5.0, 5.0]])
        types = np.array([TYPE_B, TYPE_B, TYPE_B, TYPE_A])
        n = local_B_counts(pos, types, L, r_c=1.0)
        # bead 0 sees bead 1 across the boundary (distance 0.3); bead 3 (A) sees beads 0 and 1.
        assert n.tolist() == [1, 1, 0, 2]

    def test_no_B_gives_zeros(self):
        pos = np.random.default_rng(0).uniform(0, 5, size=(20, 3))
        assert np.all(local_B_counts(pos, np.full(20, TYPE_A), 5.0, 1.5) == 0)


class TestMarkDynamics:
    def test_validation(self):
        with pytest.raises(ValueError):
            MarkDynamics(k_off=-1, k_on=0.1)
        with pytest.raises(ValueError):
            MarkDynamics(k_off=0.1, k_on=0.1, interval_steps=0)
        p = MarkDynamics(k_off=0.3, k_on=0.1)
        assert abs(p.equilibrium_fraction_without_feedback() - 0.25) < 1e-12

    def test_hill_limits(self):
        assert hill(np.array([0.0]), 6.0, 2.0)[0] == 0.0
        assert abs(hill(np.array([6.0]), 6.0, 2.0)[0] - 0.5) < 1e-12
        assert hill(np.array([1e6]), 6.0, 2.0)[0] > 0.999

    def test_no_feedback_relaxes_to_kinetic_fraction(self):
        rng = np.random.default_rng(1)
        p = MarkDynamics(k_off=0.02, k_on=0.01, k_fb=0.0)
        types = np.full(20000, TYPE_A, dtype=np.int8)
        n_B = np.zeros(types.size)
        for _ in range(2000):
            types, _, _ = step_marks(types, n_B, p, dt_marks=1.0, rng=rng)
        assert abs(np.mean(types == TYPE_B) - 1.0 / 3.0) < 0.01

    def test_feedback_raises_writing_where_B_is_dense(self):
        rng = np.random.default_rng(2)
        p = MarkDynamics(k_off=0.0, k_on=0.0, k_fb=1.0, n_half=2.0)
        types = np.full(10000, TYPE_A, dtype=np.int8)
        n_B = np.concatenate([np.zeros(5000), np.full(5000, 20.0)])
        new, n_on, n_off = step_marks(types, n_B, p, dt_marks=1.0, rng=rng)
        assert np.all(new[:5000] == TYPE_A)
        assert 0.6 < np.mean(new[5000:] == TYPE_B) < 0.66  # 1 - exp(-0.99)
        assert n_off == 0 and n_on == np.sum(new == TYPE_B)

    def test_mean_field_fixed_point_without_feedback(self):
        p = MarkDynamics(k_off=0.3, k_on=0.1)
        assert abs(mean_field_rates(p, np.array([0.25]), 6.0)[0]) < 1e-12
        assert mean_field_rates(p, np.array([0.1]), 6.0)[0] > 0 > mean_field_rates(p, np.array([0.5]), 6.0)[0]

    def test_mean_field_bistability_with_feedback(self):
        # Strong feedback and turnover create a low and a high B fixed point.
        p = MarkDynamics(k_off=0.05, k_on=0.001, k_fb=0.2, n_half=3.0, hill=3.0)
        f = np.linspace(0.0, 1.0, 2001)
        r = mean_field_rates(p, f, coordination=8.0)
        crossings = np.count_nonzero(np.diff(np.sign(r)) != 0)
        assert crossings == 3


class TestMarkRecorderAndMemory:
    def test_recorder_roundtrip(self, tmp_path):
        p = MarkDynamics(k_off=0.1, k_on=0.1, interval_steps=100)
        rec = MarkRecorder(p, 5)
        rec.observe(100, np.array([0, 1, 0, 1, 1]), 1, 0, np.zeros(5))
        rec.observe(200, np.array([0, 0, 0, 1, 1]), 1, 0, np.ones(5))
        path = rec.save(str(tmp_path / "marks.npz"), production_steps=200)
        with np.load(path) as d:
            assert d["types"].shape == (2, 5) and d["types"].dtype == np.int8
            assert np.allclose(d["f_B"], [0.4, 0.6])
            assert float(d["param_k_off"]) == 0.1 and int(d["production_steps"]) == 200

    def test_autocorrelation_and_persistence(self):
        rng = np.random.default_rng(3)
        T, N = 4000, 200
        # Telegraph process per site with switching probability p per frame.
        p = 0.05
        m = np.empty((T, N), dtype=np.int8)
        m[0] = rng.random(N) < 0.5
        for t in range(1, T):
            flip = rng.random(N) < p
            m[t] = np.where(flip, 1 - m[t - 1], m[t - 1])
        C = mark_autocorrelation(m, max_lag=40)
        assert abs(C[0] - 1.0) < 1e-9
        expected = (1 - 2 * p) ** np.arange(41)
        np.testing.assert_allclose(C, expected, atol=0.05)
        frozen = np.full((50, 10), TYPE_B, dtype=np.int8)
        frozen[:, :3] = TYPE_A
        assert abs(persistence_fraction(frozen) - 0.7) < 1e-12


class TestClusters:
    def test_pairs_and_components(self):
        L = 10.0
        pos = np.array([[0.1, 0, 0], [9.9, 0, 0], [5, 5, 5], [5.8, 5, 5], [2, 2, 2]])
        pairs = contact_pairs(pos, L, 1.0)
        assert sorted(map(tuple, pairs.tolist())) == [(0, 1), (2, 3)]
        labels = cluster_labels(5, pairs)
        assert labels[0] == labels[1] and labels[2] == labels[3] and len(set(labels)) == 3

    def test_condensation_summary_dispersed_vs_condensed(self):
        rng = np.random.default_rng(4)
        L = 20.0
        n = 400
        types = np.where(np.arange(n) < 200, TYPE_B, TYPE_A)
        dispersed = rng.uniform(0, L, size=(n, 3))
        s1 = condensation_summary(dispersed, types, L, r_c=1.0, n_dense=4)
        # Put every B bead inside a 3-sigma ball: one cluster, everything dense.
        condensed = dispersed.copy()
        condensed[:200] = 10.0 + rng.uniform(-1.5, 1.5, size=(200, 3))
        s2 = condensation_summary(condensed, types, L, r_c=1.5, n_dense=4)
        assert s1["f_B"] == 0.5 and s2["f_B"] == 0.5
        assert s2["largest_B_cluster_fraction"] > 0.95 > s1["largest_B_cluster_fraction"]
        assert s2["dense_B_fraction"] > 0.9 > s1["dense_B_fraction"]
        assert s2["n_BB_contacts"] > s1["n_BB_contacts"]

    def test_no_B_beads(self):
        s = condensation_summary(np.random.default_rng(5).uniform(0, 5, size=(30, 3)), np.full(30, TYPE_A), 5.0)
        assert s["n_B"] == 0 and np.isnan(s["largest_B_cluster_fraction"])


@needs_openmm
class TestEngineIntegration:
    def test_update_types_changes_energy_in_context(self):
        import openmm as mm

        from melt.integrator import build_openmm_system, make_context, make_langevin_integrator, update_types_in_context
        from melt.model import MeltParams

        rng = np.random.default_rng(6)
        n_chains, N = 4, 6
        L = 8.0
        pos = rng.uniform(1.0, L - 1.0, size=(n_chains * N, 3))
        types = np.full(n_chains * N, TYPE_A, dtype=np.int8)
        mp = MeltParams(n_chains=n_chains, chain_length=N, box_size=L, lj_eps_AA=0.0, lj_eps_BB=2.0, lj_eps_AB=0.0, temperature=100.0)
        system = build_openmm_system(types, mp)
        ctx = make_context(system, make_langevin_integrator(100.0, 1.0, 0.005, seed=1), pos, L, platform_name="Reference")
        e_all_A = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy()._value
        assert abs(e_all_A) < 1e-9  # no attraction among A beads
        update_types_in_context(system, ctx, np.full(n_chains * N, TYPE_B, dtype=np.int8))
        e_all_B = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy()._value
        assert e_all_B < -1e-3  # B-B attraction now active

    def test_dynamic_marks_run_writes_outputs(self, tmp_path):
        from types import SimpleNamespace

        from melt.run import execute

        cfg = SimpleNamespace(
            out=str(tmp_path), run_id="marks", seed=5, n_chains=6, chain_length=8, box_size=8.0,
            sequence="correlated", f_A=0.7, block_length=4, kappa=0.5, pi=0.99,
            bond_k=200.0, bond_r0=1.0, lj_eps_AA=0.0, lj_eps_BB=1.5, lj_eps_AB=0.0,
            lj_eps_core=1.0, lj_sigma=1.0, lj_cutoff=2.5, temperature=1.0,
            T_equilibrate=1.0, T_quench=1.0, friction=1.0, dt=0.005,
            n_steps=600, equilibration=100, snapshot_interval=200, platform="Reference",
            compute_density=False, grid_size=8, save_trajectory=True, save_density_grids=False,
            record_modes=True, mode_interval=100, mode_q_max=2.0,
            marks_dynamic=True, mark_k_off=5.0, mark_k_on=5.0, mark_k_fb=0.0, mark_interval=100,
            condensation_interval=200,
        )
        rd = execute(cfg)
        with np.load(os.path.join(rd, "marks.npz")) as d:
            assert d["steps"].tolist() == [100, 200, 300, 400, 500, 600]
            assert d["types"].shape == (6, 48)
            assert (d["n_on"].sum() + d["n_off"].sum()) > 0  # marks actually flipped at these rates
            f_B = d["f_B"]
        with np.load(os.path.join(rd, "mode_amplitudes.npz")) as d:
            assert d["n_B_t"].shape == (6,)
            assert not np.all(d["n_B_t"] == d["n_B_t"][0])
        rows = list(csv.DictReader(open(os.path.join(rd, "condensation.csv"))))
        assert [int(r["step"]) for r in rows] == [200, 400, 600]
        assert all("largest_B_cluster_fraction" in r for r in rows)
        meta = json.load(open(os.path.join(rd, "meta.json")))
        assert meta["marks"]["dynamic"] is True and meta["marks"]["k_off"] == 5.0
        # The saved trajectory carries the final marks, not the initial ones.
        with np.load(os.path.join(rd, "trajectory.npz")) as d:
            assert d["types"].shape == (48,)
        assert 0.0 < float(f_B[-1]) < 1.0
