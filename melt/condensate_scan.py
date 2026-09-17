"""Condensation campaigns with asymmetric energetics and free volume: protein and chromatin.

Both applications use "B loves B, A-B neutral" energetics (``eps_BB`` varied, ``eps_AA`` and
``eps_AB`` small) in a box with free volume (bead density well below the melt), so B-rich
material can condense and expel the implicit solvent while A stays open.

* protein mode (static marks): the sequence is quenched; the campaign scans the B-B
  attraction against the sequence correlation ``kappa`` (blockiness of the hydrophobic
  pattern) and asks when and how the chains condense.
* chromatin mode (``--marks-dynamic``): marks turn over at ``k_off`` and are written at
  ``k_on + k_fb H(n_B)``; the campaign scans turnover and feedback gain at a fixed
  ``eps_BB`` chosen near the condensation boundary, and asks whether transient blobs are
  stabilized by reader-writer feedback.

Run ids encode the varied parameters; the design is hashed into the campaign manifest.
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
from typing import Any

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
BASE_OUTPUTS = ("meta.json", "snapshots.csv", "mode_amplitudes.npz", "condensation.csv", "planned_run.json")


@dataclass(frozen=True)
class RunSpec:
    index: int
    run_id: str
    n_chains: int
    chain_length: int
    box_size: float
    bead_density: float
    eps_BB: float
    kappa: float
    f_A: float
    marks_dynamic: bool
    k_off: float
    k_on: float
    k_fb: float
    seed: int


def _num(value: float, decimals: int = 4) -> str:
    return f"{float(value):.{decimals}g}".replace("-", "m").replace(".", "p")


def make_run_id(spec_fields: dict[str, Any]) -> str:
    parts = [f"epsBB{_num(spec_fields['eps_BB'])}", f"kappa{_num(spec_fields['kappa'])}", f"fA{_num(spec_fields['f_A'])}"]
    if spec_fields["marks_dynamic"]:
        parts += [f"koff{_num(spec_fields['k_off'])}", f"kfb{_num(spec_fields['k_fb'])}"]
    parts.append(f"seed{int(spec_fields['seed'])}")
    return "_".join(parts)


def build_plan(args: argparse.Namespace) -> list[RunSpec]:
    M, N = int(args.n_chains), int(args.chain_length)
    density = float(args.bead_density)
    L = (M * N / density) ** (1.0 / 3.0)
    dynamic = bool(args.marks_dynamic)
    k_offs = [float(x) for x in args.k_offs] if dynamic else [0.0]
    k_fbs = [float(x) for x in args.k_fbs] if dynamic else [0.0]
    plan: list[RunSpec] = []
    for eps_BB in args.eps_BBs:
        for kappa in args.kappas:
            for f_A in args.f_As:
                for k_off in k_offs:
                    for k_fb in k_fbs:
                        f_B = 1.0 - float(f_A)
                        k_on = (float(args.k_on) if args.k_on is not None else k_off * f_B / max(f_A, 1e-9)) if dynamic else 0.0
                        for seed in args.seeds:
                            fields = dict(eps_BB=float(eps_BB), kappa=float(kappa), f_A=float(f_A), marks_dynamic=dynamic,
                                          k_off=k_off, k_on=k_on, k_fb=k_fb, seed=int(seed))
                            plan.append(RunSpec(index=len(plan), run_id=make_run_id(fields), n_chains=M, chain_length=N,
                                                box_size=L, bead_density=density, **fields))
    ids = [s.run_id for s in plan]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate run ids: the grid contains repeated conditions")
    return plan


def simulation_settings(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "sequence": "correlated",
        "exact_global_composition": True,
        "sequence_balance_max_attempts": 100000,
        "block_length": 4,
        "pi": float(args.pi),
        "bond_k": 200.0,
        "bond_r0": 1.0,
        "lj_eps_AA": float(args.eps_AA),
        "lj_eps_AB": float(args.eps_AB),
        "lj_eps_core": 1.0,
        "lj_sigma": 1.0,
        "lj_cutoff": 2.5,
        "T_equilibrate": float(args.T_equilibrate),
        "T_quench": float(args.T_quench),
        "temperature": float(args.T_quench),
        "friction": 1.0,
        "dt": 0.005,
        "equilibration": int(args.equilibration),
        "n_steps": int(args.n_steps),
        "snapshot_interval": int(args.snapshot_interval),
        "compute_density": False,
        "grid_size": 32,
        "save_trajectory": bool(args.save_trajectory),
        "save_density_grids": False,
        "compute_direct_structure_factor": False,
        "record_modes": True,
        "mode_interval": int(args.mode_interval),
        "mode_q_max": float(args.mode_q_max),
        "mode_dtype": "complex64",
        "condensation_interval": int(args.condensation_interval),
        "condensation_r_c": float(args.contact_radius),
        "condensation_n_dense": int(args.n_dense),
        "mark_r_c": float(args.contact_radius),
        "mark_n_half": float(args.mark_n_half),
        "mark_hill": float(args.mark_hill),
        "mark_interval": int(args.mark_interval),
    }


def campaign_design(args: argparse.Namespace, plan: list[RunSpec]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "study": "chromatin epigenetic memory" if args.marks_dynamic else "protein sequence-programmed condensation",
        "campaign_id": args.campaign_id,
        "design": {
            "mode": "chromatin" if args.marks_dynamic else "protein",
            "eps_BBs": [float(x) for x in args.eps_BBs],
            "kappas": [float(x) for x in args.kappas],
            "f_As": [float(x) for x in args.f_As],
            "k_offs": [float(x) for x in args.k_offs] if args.marks_dynamic else [],
            "k_fbs": [float(x) for x in args.k_fbs] if args.marks_dynamic else [],
            "k_on": ("auto: k_off * f_B/f_A (steady B fraction equals the initial one without feedback)"
                     if args.k_on is None else float(args.k_on)),
            "seeds": [int(x) for x in args.seeds],
            "n_chains": int(args.n_chains),
            "chain_length": int(args.chain_length),
            "bead_density": float(args.bead_density),
            "energetics": "B-B attraction eps_BB varied; eps_AA and eps_AB fixed small (B loves B, A-B neutral); shared WCA core",
            "n_runs": len(plan),
        },
        "simulation_settings": simulation_settings(args),
        "runs": [asdict(spec) for spec in plan],
    }


def required_outputs(dynamic: bool, save_trajectory: bool) -> tuple[str, ...]:
    outputs = BASE_OUTPUTS + (("marks.npz",) if dynamic else ())
    return outputs + (("trajectory.npz",) if save_trajectory else ())


def run_config(spec: RunSpec, settings: dict[str, Any], platform: str | None) -> dict[str, Any]:
    return {
        **settings,
        "platform": platform,
        "seed": spec.seed,
        "n_chains": spec.n_chains,
        "chain_length": spec.chain_length,
        "box_size": spec.box_size,
        "kappa": spec.kappa,
        "f_A": spec.f_A,
        "lj_eps_BB": spec.eps_BB,
        "marks_dynamic": spec.marks_dynamic,
        "mark_k_off": spec.k_off,
        "mark_k_on": spec.k_on,
        "mark_k_fb": spec.k_fb,
        "skip_if_cached": False,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Condensation campaigns: protein (static marks) or chromatin (dynamic marks)")
    mode = parser.add_mutually_exclusive_group(required=True)
    for flag, help_text in (("--dry-run", "print the plan"), ("--manifest-only", "write/validate the manifest"),
                            ("--run", "execute runs"), ("--status", "report completion"), ("--analyze", "aggregate")):
        mode.add_argument(flag, action="store_true", help=help_text)
    parser.add_argument("--out", default="output/melt/condensate")
    parser.add_argument("--campaign-id", required=False, default="condensate")
    parser.add_argument("--marks-dynamic", action="store_true", help="chromatin mode: marks turn over and are written")
    parser.add_argument("--eps-BBs", dest="eps_BBs", type=float, nargs="+", default=[0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])
    parser.add_argument("--kappas", type=float, nargs="+", default=[0.0, 0.5, 1.0])
    parser.add_argument("--f-As", dest="f_As", type=float, nargs="+", default=[0.7], help="A (unmarked) fraction")
    parser.add_argument("--k-offs", dest="k_offs", type=float, nargs="+", default=[1e-3])
    parser.add_argument("--k-fbs", dest="k_fbs", type=float, nargs="+", default=[0.0])
    parser.add_argument("--k-on", dest="k_on", type=float, default=None, help="basal writing rate (default: auto)")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4])
    parser.add_argument("--n-chains", type=int, default=144)
    parser.add_argument("--chain-length", type=int, default=40)
    parser.add_argument("--bead-density", type=float, default=0.2)
    parser.add_argument("--eps-AA", dest="eps_AA", type=float, default=0.0)
    parser.add_argument("--eps-AB", dest="eps_AB", type=float, default=0.0)
    parser.add_argument("--pi", type=float, default=0.99)
    parser.add_argument("--T-equilibrate", dest="T_equilibrate", type=float, default=3.0)
    parser.add_argument("--T-quench", dest="T_quench", type=float, default=1.0)
    parser.add_argument("--equilibration", type=int, default=30000)
    parser.add_argument("--n-steps", type=int, default=1000000)
    parser.add_argument("--snapshot-interval", type=int, default=2000)
    parser.add_argument("--mode-interval", type=int, default=400)
    parser.add_argument("--mode-q-max", type=float, default=1.5)
    parser.add_argument("--condensation-interval", type=int, default=2000)
    parser.add_argument("--contact-radius", type=float, default=1.5)
    parser.add_argument("--n-dense", type=int, default=6)
    parser.add_argument("--mark-n-half", type=float, default=6.0)
    parser.add_argument("--mark-hill", type=float, default=2.0)
    parser.add_argument("--mark-interval", type=int, default=200)
    parser.add_argument("--save-trajectory", action="store_true")
    parser.add_argument("--platform", default=None)
    parser.add_argument("--run-index", type=int, action="append", default=None)
    parser.add_argument("--shard-index", type=int, default=None)
    parser.add_argument("--shard-count", type=int, default=None)
    parser.add_argument("--parallel", type=int, default=1)
    parser.add_argument("--threads-per-run", type=int, default=None)
    parser.add_argument("--cuda-devices", type=str, default=None)
    parser.add_argument("--recover-interrupted", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument("--allow-incomplete-analysis", action="store_true")
    parser.add_argument("--window-fraction", type=float, default=0.5)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    n = int(args.n_steps)
    for name in ("snapshot_interval", "mode_interval", "condensation_interval", "mark_interval"):
        v = int(getattr(args, name))
        if v < 1 or n % v:
            raise ValueError(f"{name} must be positive and divide n_steps")
    if int(args.condensation_interval) % int(args.snapshot_interval):
        raise ValueError("condensation_interval must be a multiple of snapshot_interval")
    if float(args.bead_density) <= 0 or float(args.bead_density) > 1.0:
        raise ValueError("bead_density must lie in (0,1]")
    if any(not 0.0 < float(f) < 1.0 for f in args.f_As):
        raise ValueError("f_As must lie strictly inside (0,1)")
    n_beads = int(args.n_chains) * int(args.chain_length)
    for f in args.f_As:
        if abs(float(f) * n_beads - round(float(f) * n_beads)) > 1e-6:
            raise ValueError(f"f_A={f} times {n_beads} beads is not an integer; exact composition needs one")
    if any(float(k) < 0 for k in args.eps_BBs) or any(not 0 <= float(k) <= 1 for k in args.kappas):
        raise ValueError("eps_BBs must be nonnegative and kappas in [0,1]")


def _child_command(args: argparse.Namespace, index: int) -> list[str]:
    cmd = [sys.executable, "-m", "melt.condensate_scan", "--run", "--out", args.out, "--campaign-id", args.campaign_id,
           "--run-index", str(index)]
    if args.marks_dynamic:
        cmd.append("--marks-dynamic")
    for flag, values in (("--eps-BBs", args.eps_BBs), ("--kappas", args.kappas), ("--f-As", args.f_As),
                         ("--k-offs", args.k_offs), ("--k-fbs", args.k_fbs), ("--seeds", args.seeds)):
        cmd += [flag, *[repr(float(v)) if flag != "--seeds" else str(int(v)) for v in values]]
    if args.k_on is not None:
        cmd += ["--k-on", repr(float(args.k_on))]
    for flag, attr in (("--n-chains", "n_chains"), ("--chain-length", "chain_length"), ("--bead-density", "bead_density"),
                       ("--eps-AA", "eps_AA"), ("--eps-AB", "eps_AB"), ("--pi", "pi"), ("--T-equilibrate", "T_equilibrate"),
                       ("--T-quench", "T_quench"), ("--equilibration", "equilibration"), ("--n-steps", "n_steps"),
                       ("--snapshot-interval", "snapshot_interval"), ("--mode-interval", "mode_interval"),
                       ("--mode-q-max", "mode_q_max"), ("--condensation-interval", "condensation_interval"),
                       ("--contact-radius", "contact_radius"), ("--n-dense", "n_dense"), ("--mark-n-half", "mark_n_half"),
                       ("--mark-hill", "mark_hill"), ("--mark-interval", "mark_interval")):
        cmd += [flag, repr(getattr(args, attr))]
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


def status_report(campaign_dir: Path, plan: list[RunSpec], outputs: tuple[str, ...], commit: str | None) -> dict[str, Any]:
    missing = [s.index for s in plan if not completion_valid(campaign_dir / "runs" / s.run_id, asdict(s), outputs, commit)]
    return {"n_runs": len(plan), "n_completed": len(plan) - len(missing), "n_missing": len(missing), "missing_indices": missing}


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    validate_args(args)
    repo_root = Path(__file__).resolve().parents[1]
    plan = build_plan(args)
    design = campaign_design(args, plan)
    digest = design_hash(design, DESIGN_KEYS)
    outputs = required_outputs(bool(args.marks_dynamic), bool(args.save_trajectory))
    selected = select_plan(plan, args.run_index)
    indices = set(select_shard([s.index for s in selected], args.shard_index, args.shard_count))
    selected = [s for s in selected if s.index in indices]

    if args.dry_run:
        print(json.dumps({"design_sha256": digest, "n_runs": len(plan), "box_size": plan[0].box_size if plan else None,
                          "runs": [asdict(s) for s in plan]}, indent=2))
        return 0
    campaign_dir = Path(args.out).resolve() / args.campaign_id
    if args.status or args.analyze:
        with (campaign_dir / "manifest.json").open() as handle:
            manifest = json.load(handle)
        validate_manifest(manifest, DESIGN_KEYS)
        commit = campaign_commit(manifest)
        if args.status:
            print(json.dumps(status_report(campaign_dir, plan, outputs, commit), indent=2))
            return 0
        from .condensate_analysis import aggregate_campaign

        result = aggregate_campaign(campaign_dir, require_complete=not args.allow_incomplete_analysis,
                                    window_fraction=float(args.window_fraction))
        print(json.dumps({k: str(v) for k, v in result.items()}, indent=2, sort_keys=True))
        return 0

    manifest = ensure_campaign(campaign_dir, design, DESIGN_KEYS, repo_root, require_clean_git=not args.allow_dirty)
    commit = campaign_commit(manifest)
    print(f"campaign={args.campaign_id} design={manifest['design_sha256']} runs={len(plan)} L={plan[0].box_size:.3f} directory={campaign_dir}")
    if args.manifest_only:
        return 0
    print(f"selected_runs={len(selected)}")
    if int(args.parallel) > 1 and len(selected) > 1:
        pending = [s.index for s in selected if not completion_valid(campaign_dir / "runs" / s.run_id, asdict(s), outputs, commit)]
        print(f"cached={len(selected) - len(pending)} pending={len(pending)} parallel={args.parallel}")
        failures = run_parallel(args, pending)
        print(json.dumps(status_report(campaign_dir, plan, outputs, commit), indent=2))
        return 1 if failures else 0
    settings = simulation_settings(args)
    for position, spec in enumerate(selected, 1):
        print(f"[{position}/{len(selected)}] index={spec.index} {spec.run_id}", flush=True)
        status, output = run_staged(spec.run_id, asdict(spec), run_config(spec, settings, args.platform), outputs,
                                    campaign_dir, repo_root, commit, args.recover_interrupted)
        print(f"  {status}: {output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
