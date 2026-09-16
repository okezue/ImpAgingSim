"""Restartable ``eps_AB`` sweep at fixed sequence correlation: locate the mixed -> demixed transition.

Everything except the A-B attraction ``eps_AB`` is held at the production values
(144 x 40 beads, L = 22 sigma, kappa = 0.5, pi = 0.99, f_A = 1/2, T* = 0.7).  With
``eps_AA = eps_BB = 1`` the incompatibility is ``1 - eps_AB``: ``eps_AB = 1`` is the
chi = 0 reference where every pair interacts identically and the melt must mix, and the
grid descends to the strongly demixed production value ``eps_AB = 0.1``.

Every run records exact box-mode amplitudes ``rho_A(q,t)``, ``rho_B(q,t)`` so both the
static fluctuation spectrum and the dynamic structure factor come out of the same data.
Run plans are hashed into a manifest; runs are staged, hashed and promoted atomically, so
the sweep can be sharded over a scheduler array or many local processes and re-invoked.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

import numpy as np

from .campaign import (
    campaign_commit,
    completion_valid,
    design_hash,
    ensure_campaign,
    run_staged,
    select_shard,
    validate_manifest,
)
from .fixed_density_size_scan import select_plan

DESIGN_KEYS = ("schema_version", "study", "campaign_id", "design", "simulation_settings", "runs")
DEFAULT_EPS_ABS = tuple(float(f"{x:.3f}") for x in np.arange(1.0, 0.1 - 1e-9, -0.025))
DEFAULT_SEEDS = (1, 2, 3, 4)
DEFAULT_SIZES = (144,)
REQUIRED_OUTPUTS = ("meta.json", "snapshots.csv", "structure_factor.npz", "mode_amplitudes.npz", "planned_run.json")


@dataclass(frozen=True)
class RunSpec:
    index: int
    run_id: str
    n_chains: int
    chain_length: int
    n_beads: int
    box_size: float
    bead_density: float
    grid_size: int
    kappa: float
    pi: float
    eps_AB: float
    seed: int


def _num_text(value: float, decimals: int) -> str:
    return f"{float(value):.{decimals}f}".rstrip("0").rstrip(".").replace(".", "p")


def make_run_id(n_chains: int, chain_length: int, kappa: float, pi: float, eps_AB: float, seed: int) -> str:
    pi_text = f"{float(pi):.3f}".replace(".", "p")
    return (
        f"M{int(n_chains):04d}_N{int(chain_length):02d}_kappa{_num_text(kappa, 3)}_"
        f"pi{pi_text}_epsAB{_num_text(eps_AB, 4)}_seed{int(seed)}"
    )


def build_plan(
    eps_ABs: Iterable[float] = DEFAULT_EPS_ABS,
    seeds: Iterable[int] = DEFAULT_SEEDS,
    sizes: Iterable[int] = DEFAULT_SIZES,
    kappa: float = 0.5,
    pi: float = 0.99,
    chain_length: int = 40,
    reference_chains: int = 144,
    reference_box_size: float = 22.0,
    reference_grid_size: int = 56,
) -> list[RunSpec]:
    eps_t = tuple(float(x) for x in eps_ABs)
    seeds_t = tuple(int(x) for x in seeds)
    sizes_t = tuple(int(x) for x in sizes)
    if not eps_t or any(x < 0.0 for x in eps_t) or len(set(eps_t)) != len(eps_t):
        raise ValueError("eps_ABs must be nonempty, nonnegative and unique")
    if not seeds_t or len(set(seeds_t)) != len(seeds_t):
        raise ValueError("seeds must be nonempty and unique")
    if not sizes_t or any(x < 2 or x % 2 for x in sizes_t) or len(set(sizes_t)) != len(sizes_t):
        raise ValueError("sizes must be nonempty, even, positive and unique")
    if not 0.0 <= kappa <= 1.0 or not 0.0 <= pi <= 1.0:
        raise ValueError("kappa and pi must lie in [0,1]")
    plan: list[RunSpec] = []
    for M in sizes_t:
        scale = (M / float(reference_chains)) ** (1.0 / 3.0)
        L = float(reference_box_size) * scale
        G = max(2, int(round(reference_grid_size * scale)))
        n_beads = M * int(chain_length)
        for eps in eps_t:
            for seed in seeds_t:
                plan.append(
                    RunSpec(
                        index=len(plan),
                        run_id=make_run_id(M, chain_length, kappa, pi, eps, seed),
                        n_chains=M,
                        chain_length=int(chain_length),
                        n_beads=n_beads,
                        box_size=L,
                        bead_density=n_beads / L**3,
                        grid_size=G,
                        kappa=float(kappa),
                        pi=float(pi),
                        eps_AB=eps,
                        seed=seed,
                    )
                )
    return plan


def simulation_settings(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "sequence": "correlated",
        "exact_global_composition": True,
        "sequence_balance_max_attempts": int(args.sequence_balance_max_attempts),
        "f_A": 0.5,
        "block_length": 4,
        "bond_k": float(args.bond_k),
        "bond_r0": 1.0,
        "lj_eps_AA": float(args.lj_eps_like),
        "lj_eps_BB": float(args.lj_eps_like),
        "lj_eps_core": 1.0,
        "lj_sigma": 1.0,
        "lj_cutoff": 2.5,
        "T_equilibrate": float(args.T_equilibrate),
        "T_quench": float(args.T_quench),
        "temperature": float(args.T_quench),
        "friction": 1.0,
        "dt": float(args.dt),
        "equilibration": int(args.equilibration),
        "n_steps": int(args.n_steps),
        "snapshot_interval": int(args.snapshot_interval),
        "compute_density": True,
        "save_trajectory": bool(args.save_trajectory),
        "save_density_grids": False,
        "compute_direct_structure_factor": False,
        "record_modes": True,
        "mode_interval": int(args.mode_interval),
        "mode_q_max": float(args.mode_q_max),
        "mode_dtype": "complex64",
    }


def campaign_design(args: argparse.Namespace, plan: list[RunSpec]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "study": "eps_AB incompatibility sweep at fixed sequence correlation",
        "campaign_id": args.campaign_id,
        "design": {
            "eps_ABs": [float(x) for x in args.eps_ABs],
            "seeds": [int(x) for x in args.seeds],
            "sizes": [int(x) for x in args.sizes],
            "kappa": float(args.kappa),
            "pi": float(args.pi),
            "chain_length": int(args.chain_length),
            "reference_chains": 144,
            "reference_box_size": 22.0,
            "box_scaling": "L(M)=22*(M/144)^(1/3)",
            "reference_grid_size": 56,
            "incompatibility_axis": "eps_like - eps_AB with eps_AA=eps_BB=eps_like; eps_AB=eps_like is the chi=0 reference",
            "n_runs": len(plan),
        },
        "simulation_settings": simulation_settings(args),
        "runs": [asdict(spec) for spec in plan],
    }


def required_outputs(save_trajectory: bool) -> tuple[str, ...]:
    return (*REQUIRED_OUTPUTS, "trajectory.npz") if save_trajectory else REQUIRED_OUTPUTS


def run_config(spec: RunSpec, settings: dict[str, Any], platform: str | None) -> dict[str, Any]:
    return {
        **settings,
        "platform": platform,
        "seed": spec.seed,
        "n_chains": spec.n_chains,
        "chain_length": spec.chain_length,
        "box_size": spec.box_size,
        "grid_size": spec.grid_size,
        "kappa": spec.kappa,
        "pi": spec.pi,
        "lj_eps_AB": spec.eps_AB,
        "skip_if_cached": False,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Restartable eps_AB sweep at fixed kappa")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="print the plan; write nothing")
    mode.add_argument("--manifest-only", action="store_true", help="write/validate the manifest; run nothing")
    mode.add_argument("--run", action="store_true", help="execute selected runs")
    mode.add_argument("--status", action="store_true", help="report completed/missing runs")
    mode.add_argument("--analyze", action="store_true", help="aggregate completed runs")
    parser.add_argument("--out", default="output/melt/epsab")
    parser.add_argument("--campaign-id", default="epsab_kappa05")
    parser.add_argument("--eps-ABs", dest="eps_ABs", type=float, nargs="+", default=list(DEFAULT_EPS_ABS))
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    parser.add_argument("--sizes", type=int, nargs="+", default=list(DEFAULT_SIZES))
    parser.add_argument("--kappa", type=float, default=0.5)
    parser.add_argument("--pi", type=float, default=0.99)
    parser.add_argument("--chain-length", type=int, default=40)
    parser.add_argument("--lj-eps-like", dest="lj_eps_like", type=float, default=1.0)
    parser.add_argument("--bond-k", type=float, default=200.0)
    parser.add_argument("--T-equilibrate", dest="T_equilibrate", type=float, default=5.0)
    parser.add_argument("--T-quench", dest="T_quench", type=float, default=0.7)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--equilibration", type=int, default=30000)
    parser.add_argument("--n-steps", type=int, default=250000)
    parser.add_argument("--snapshot-interval", type=int, default=2000)
    parser.add_argument("--mode-interval", type=int, default=200)
    parser.add_argument("--mode-q-max", type=float, default=1.5)
    parser.add_argument("--sequence-balance-max-attempts", type=int, default=100000)
    parser.add_argument("--save-trajectory", action="store_true")
    parser.add_argument("--platform", default=None, help="OpenMM platform: CUDA, OpenCL, CPU or Reference")
    parser.add_argument("--run-index", type=int, action="append", default=None)
    parser.add_argument("--shard-index", type=int, default=None)
    parser.add_argument("--shard-count", type=int, default=None)
    parser.add_argument("--parallel", type=int, default=1, help="run this many runs concurrently as subprocesses")
    parser.add_argument("--threads-per-run", type=int, default=None, help="OPENMM_CPU_THREADS for each subprocess")
    parser.add_argument("--cuda-devices", type=str, default=None,
                        help="comma-separated GPU ids to round-robin over parallel subprocesses, e.g. 0,1,2,3")
    parser.add_argument("--recover-interrupted", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true", help="skip the clean-git-tree guard (smoke tests only)")
    parser.add_argument("--allow-incomplete-analysis", action="store_true")
    parser.add_argument("--window-fraction", type=float, default=0.5, help="trailing fraction of production analyzed")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if int(args.n_steps) < 1 or int(args.snapshot_interval) < 1 or int(args.mode_interval) < 1:
        raise ValueError("n_steps, snapshot_interval and mode_interval must be positive")
    if int(args.n_steps) % int(args.snapshot_interval) or int(args.n_steps) % int(args.mode_interval):
        raise ValueError("snapshot_interval and mode_interval must divide n_steps")
    if int(args.equilibration) < 0:
        raise ValueError("equilibration must be nonnegative")
    if float(args.mode_q_max) <= 0.0:
        raise ValueError("mode_q_max must be positive")
    if int(args.parallel) < 1:
        raise ValueError("parallel must be positive")
    if not 0.0 < float(args.window_fraction) <= 1.0:
        raise ValueError("window_fraction must lie in (0,1]")


def _spec_dict(spec: RunSpec) -> dict[str, Any]:
    return asdict(spec)


def status_report(campaign_dir: Path, plan: list[RunSpec], save_trajectory: bool, commit: str | None) -> dict[str, Any]:
    done, missing = [], []
    for spec in plan:
        if completion_valid(campaign_dir / "runs" / spec.run_id, _spec_dict(spec), required_outputs(save_trajectory), commit):
            done.append(spec.index)
        else:
            missing.append(spec.index)
    return {"n_runs": len(plan), "n_completed": len(done), "n_missing": len(missing), "missing_indices": missing}


def _child_command(args: argparse.Namespace, index: int) -> list[str]:
    cmd = [sys.executable, "-m", "melt.epsab_scan", "--run", "--out", args.out, "--campaign-id", args.campaign_id,
           "--run-index", str(index), "--eps-ABs", *[repr(float(x)) for x in args.eps_ABs],
           "--seeds", *[str(int(x)) for x in args.seeds], "--sizes", *[str(int(x)) for x in args.sizes],
           "--kappa", repr(float(args.kappa)), "--pi", repr(float(args.pi)), "--chain-length", str(args.chain_length),
           "--lj-eps-like", repr(float(args.lj_eps_like)), "--bond-k", repr(float(args.bond_k)),
           "--T-equilibrate", repr(float(args.T_equilibrate)), "--T-quench", repr(float(args.T_quench)),
           "--dt", repr(float(args.dt)), "--equilibration", str(args.equilibration), "--n-steps", str(args.n_steps),
           "--snapshot-interval", str(args.snapshot_interval), "--mode-interval", str(args.mode_interval),
           "--mode-q-max", repr(float(args.mode_q_max)),
           "--sequence-balance-max-attempts", str(args.sequence_balance_max_attempts)]
    if args.save_trajectory:
        cmd.append("--save-trajectory")
    if args.platform:
        cmd += ["--platform", args.platform]
    if args.recover_interrupted:
        cmd.append("--recover-interrupted")
    if args.allow_dirty:
        cmd.append("--allow-dirty")
    return cmd


def run_parallel(args: argparse.Namespace, indices: list[int]) -> int:
    env = dict(os.environ)
    if args.threads_per_run:
        env["OPENMM_CPU_THREADS"] = str(int(args.threads_per_run))
        env.setdefault("OMP_NUM_THREADS", "1")
        env.setdefault("OPENBLAS_NUM_THREADS", "1")
    devices = [d.strip() for d in (args.cuda_devices or "").split(",") if d.strip()]
    failures = 0

    def worker(slot_index: tuple[int, int]) -> tuple[int, int, str]:
        slot, index = slot_index
        child_env = dict(env)
        if devices:
            # Round-robin the workers over the listed GPUs; each child sees exactly one device.
            child_env["CUDA_VISIBLE_DEVICES"] = devices[slot % len(devices)]
        proc = subprocess.run(_child_command(args, index), env=child_env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return index, proc.returncode, proc.stdout

    with ThreadPoolExecutor(max_workers=int(args.parallel)) as pool:
        for index, code, output in pool.map(worker, list(enumerate(indices))):
            tail = output.strip().splitlines()[-1] if output.strip() else ""
            print(f"[index={index}] exit={code} {tail}", flush=True)
            if code != 0:
                failures += 1
                print(output, file=sys.stderr, flush=True)
    return failures


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    validate_args(args)
    repo_root = Path(__file__).resolve().parents[1]
    plan = build_plan(args.eps_ABs, args.seeds, args.sizes, args.kappa, args.pi, args.chain_length)
    design = campaign_design(args, plan)
    digest = design_hash(design, DESIGN_KEYS)
    selected = select_plan(plan, args.run_index)
    indices = select_shard([spec.index for spec in selected], args.shard_index, args.shard_count)
    selected = [spec for spec in selected if spec.index in set(indices)]

    if args.dry_run:
        print(json.dumps({"design_sha256": digest, "n_runs": len(plan), "runs": [asdict(s) for s in plan]}, indent=2))
        return 0

    campaign_dir = Path(args.out).resolve() / args.campaign_id
    if args.analyze or args.status:
        manifest_path = campaign_dir / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"campaign manifest does not exist: {manifest_path}")
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        validate_manifest(manifest, DESIGN_KEYS)
        commit = campaign_commit(manifest)
        if args.status:
            print(json.dumps(status_report(campaign_dir, plan, args.save_trajectory, commit), indent=2))
            return 0
        from .epsab_analysis import aggregate_campaign

        outputs = aggregate_campaign(
            campaign_dir,
            require_complete=not args.allow_incomplete_analysis,
            window_fraction=float(args.window_fraction),
        )
        print(json.dumps({key: str(value) for key, value in outputs.items()}, indent=2, sort_keys=True))
        return 0

    manifest = ensure_campaign(campaign_dir, design, DESIGN_KEYS, repo_root, require_clean_git=not args.allow_dirty)
    commit = campaign_commit(manifest)
    print(f"campaign={args.campaign_id} design={manifest['design_sha256']} runs={len(plan)} directory={campaign_dir}")
    if args.manifest_only:
        print(json.dumps({"design_sha256": digest, "n_runs": len(plan), "runs": [asdict(s) for s in plan]}, indent=2))
        return 0

    print(f"selected_runs={len(selected)}")
    if int(args.parallel) > 1 and len(selected) > 1:
        pending = [s.index for s in selected
                   if not completion_valid(campaign_dir / "runs" / s.run_id, _spec_dict(s), required_outputs(args.save_trajectory), commit)]
        print(f"cached={len(selected)-len(pending)} pending={len(pending)} parallel={args.parallel}")
        failures = run_parallel(args, pending)
        print(json.dumps(status_report(campaign_dir, plan, args.save_trajectory, commit), indent=2))
        return 1 if failures else 0

    settings = simulation_settings(args)
    for position, spec in enumerate(selected, 1):
        print(f"[{position}/{len(selected)}] index={spec.index} {spec.run_id} eps_AB={spec.eps_AB:g} L={spec.box_size:.6f}", flush=True)
        status, output = run_staged(
            spec.run_id,
            _spec_dict(spec),
            run_config(spec, settings, args.platform),
            required_outputs(args.save_trajectory),
            campaign_dir,
            repo_root,
            commit,
            args.recover_interrupted,
        )
        print(f"  {status}: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
