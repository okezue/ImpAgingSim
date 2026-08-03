from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import socket
import subprocess
import sys
import time
from types import SimpleNamespace
from typing import Any, Iterable

import numpy as np


SCHEMA_VERSION = 1
CAMPAIGN_DESIGN_KEYS = (
    "schema_version",
    "study",
    "campaign_id",
    "design",
    "simulation_settings",
    "direct_structure_factor",
    "runs",
)
DEFAULT_SIZES = (144, 288, 576)
DEFAULT_KAPPAS = (0.0, 1.0)
DEFAULT_SEEDS = (1, 2, 3, 4, 5)


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
    grid_spacing: float
    kappa: float
    pi: float
    seed: int


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest_design_hash(manifest: dict[str, Any]) -> str:
    if not isinstance(manifest, dict):
        raise ValueError("campaign manifest root must be an object")
    missing = [key for key in CAMPAIGN_DESIGN_KEYS if key not in manifest]
    if missing:
        raise ValueError(f"campaign manifest lacks design fields {missing}")
    payload = {key: manifest[key] for key in CAMPAIGN_DESIGN_KEYS}
    actual = sha256_bytes(canonical_json(payload))
    recorded = manifest.get("design_sha256")
    if not isinstance(recorded, str) or recorded != actual:
        raise ValueError("campaign manifest failed its design hash")
    return actual


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        handle.write(json.dumps(value, indent=2, sort_keys=True).encode("utf-8"))
        handle.write(b"\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _git_output(repo_root: Path, args: list[str], binary: bool = False) -> bytes | str | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=repo_root,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return result.stdout if binary else result.stdout.decode("utf-8", errors="replace").strip()


def collect_provenance(repo_root: Path) -> dict[str, Any]:
    commit = _git_output(repo_root, ["rev-parse", "HEAD"])
    status = _git_output(repo_root, ["status", "--porcelain=v1"])
    diff = _git_output(repo_root, ["diff", "--binary", "HEAD"], binary=True)
    try:
        openmm_version = importlib.metadata.version("openmm")
    except importlib.metadata.PackageNotFoundError:
        openmm_version = None
    return {
        "git_commit": commit,
        "git_dirty": bool(status),
        "git_diff_sha256": sha256_bytes(diff) if diff else None,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "openmm_version": openmm_version,
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "python_executable": sys.executable,
        "argv": list(sys.argv),
    }


def _validate_clean_git_provenance(provenance: dict[str, Any], label: str) -> str:
    if not isinstance(provenance, dict):
        raise RuntimeError(f"{label} provenance is missing")
    commit = provenance.get("git_commit")
    if not isinstance(commit, str) or len(commit) != 40:
        raise RuntimeError(f"{label} provenance lacks a full git commit")
    if provenance.get("git_dirty") is not False or provenance.get("git_diff_sha256") is not None:
        raise RuntimeError(f"{label} git tree is dirty")
    return commit


def validate_creation_provenance(
    manifest: dict[str, Any], current_provenance: dict[str, Any] | None = None
) -> str:
    creation = manifest.get("creation_provenance")
    recorded_hash = manifest.get("creation_provenance_sha256")
    if not isinstance(creation, dict) or recorded_hash != sha256_bytes(canonical_json(creation)):
        raise RuntimeError("campaign creation provenance failed its hash")
    creation_commit = _validate_clean_git_provenance(creation, "campaign creation")
    if current_provenance is not None:
        current_commit = _validate_clean_git_provenance(current_provenance, "current")
        if current_commit != creation_commit:
            raise RuntimeError(
                f"current commit {current_commit} differs from campaign commit {creation_commit}"
            )
    return creation_commit


def make_run_id(n_chains: int, chain_length: int, kappa: float, pi: float, seed: int) -> str:
    k_text = f"{float(kappa):.3f}".rstrip("0").rstrip(".").replace(".", "p")
    pi_text = f"{float(pi):.3f}".replace(".", "p")
    return (
        f"M{int(n_chains):04d}_N{int(chain_length):02d}_"
        f"kappa{k_text}_pi{pi_text}_seed{int(seed)}"
    )


def build_plan(
    sizes: Iterable[int] = DEFAULT_SIZES,
    kappas: Iterable[float] = DEFAULT_KAPPAS,
    seeds: Iterable[int] = DEFAULT_SEEDS,
    chain_length: int = 40,
    reference_chains: int = 144,
    reference_box_size: float = 22.0,
    reference_grid_size: int = 56,
    pi: float = 0.99,
) -> list[RunSpec]:
    sizes_t = tuple(int(x) for x in sizes)
    kappas_t = tuple(float(x) for x in kappas)
    seeds_t = tuple(int(x) for x in seeds)
    if not sizes_t or any(x < 2 or x % 2 for x in sizes_t) or len(set(sizes_t)) != len(sizes_t):
        raise ValueError("sizes must be nonempty, positive even integers, and unique")
    if not kappas_t or any(x < 0.0 or x > 1.0 for x in kappas_t):
        raise ValueError("kappas must lie in [0,1]")
    if len(set(kappas_t)) != len(kappas_t):
        raise ValueError("kappas must be unique")
    if not seeds_t or len(set(seeds_t)) != len(seeds_t):
        raise ValueError("seeds must be nonempty and unique")
    if chain_length < 1 or reference_chains < 1 or reference_box_size <= 0.0:
        raise ValueError("chain length and reference system dimensions must be positive")
    if reference_grid_size < 2:
        raise ValueError("reference_grid_size must be at least 2")
    if not 0.0 <= pi <= 1.0:
        raise ValueError("pi must lie in [0,1]")

    plan: list[RunSpec] = []
    for M in sizes_t:
        scale = (M / float(reference_chains)) ** (1.0 / 3.0)
        L = float(reference_box_size) * scale
        G = max(2, int(round(reference_grid_size * scale)))
        n_beads = M * int(chain_length)
        rho = n_beads / L**3
        for kappa in kappas_t:
            for seed in seeds_t:
                plan.append(
                    RunSpec(
                        index=len(plan),
                        run_id=make_run_id(M, chain_length, kappa, pi, seed),
                        n_chains=M,
                        chain_length=int(chain_length),
                        n_beads=n_beads,
                        box_size=L,
                        bead_density=rho,
                        grid_size=G,
                        grid_spacing=L / G,
                        kappa=kappa,
                        pi=float(pi),
                        seed=seed,
                    )
                )
    return plan


def simulation_settings(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "sequence": "correlated",
        "exact_global_composition": True,
        "rng_stream_scheme": "numpy SeedSequence(seed).spawn(2): sequence and placement",
        "sequence_balance_max_attempts": int(args.sequence_balance_max_attempts),
        "f_A": float(args.f_A),
        "block_length": int(args.block_length),
        "bond_k": float(args.bond_k),
        "bond_r0": float(args.bond_r0),
        "lj_eps_AA": float(args.lj_eps_AA),
        "lj_eps_BB": float(args.lj_eps_BB),
        "lj_eps_AB": float(args.lj_eps_AB),
        "lj_eps_core": float(args.lj_eps_core),
        "lj_sigma": float(args.lj_sigma),
        "lj_cutoff": float(args.lj_cutoff),
        "T_equilibrate": float(args.T_equilibrate),
        "T_quench": float(args.T_quench),
        "temperature": float(args.T_quench),
        "friction": float(args.friction),
        "dt": float(args.dt),
        "equilibration": int(args.equilibration),
        "n_steps": int(args.n_steps),
        "snapshot_interval": int(args.snapshot_interval),
        "platform": args.platform,
        "compute_density": True,
        "save_trajectory": bool(args.save_trajectory),
        "save_density_grids": False,
        "compute_direct_structure_factor": True,
        "direct_q_max": float(args.direct_q_max),
        "direct_final_frames": int(args.direct_final_frames),
        "direct_early_frames": int(args.direct_early_frames),
        "direct_chunk_size": int(args.direct_chunk_size),
    }


def campaign_design(args: argparse.Namespace, plan: list[RunSpec]) -> dict[str, Any]:
    settings = simulation_settings(args)
    return {
        "schema_version": SCHEMA_VERSION,
        "study": "fixed-density finite-size comparison",
        "campaign_id": args.campaign_id,
        "design": {
            "sizes": [int(x) for x in args.sizes],
            "kappas": [float(x) for x in args.kappas],
            "seeds": [int(x) for x in args.seeds],
            "chain_length": int(args.chain_length),
            "reference_chains": int(args.reference_chains),
            "reference_box_size": float(args.reference_box_size),
            "box_scaling": "L(M)=reference_box_size*(M/reference_chains)^(1/3)",
            "reference_grid_size": int(args.reference_grid_size),
            "grid_scaling": "G(M)=round(reference_grid_size*(M/reference_chains)^(1/3))",
            "pi": float(args.pi),
            "sequence_ensemble": (
                "independently drawn full M-chain sets conditioned by rejection sampling on "
                "exactly M*N/2 A beads; chains need not be unique and no complement pairs are constructed"
            ),
            "exact_global_f_A": 0.5,
            "n_runs": len(plan),
        },
        "simulation_settings": settings,
        "direct_structure_factor": {
            "mode_definition": "q=2*pi*h/L, h in Z^3, 0<|q|<=direct_q_max",
            "partial_normalization": "S_ab(q)=Re[rho_a(q)rho_b(q)*]/N_total",
            "composition_channel": "S_psi_psi^(N)=S_CC=S_AA+S_BB-2*S_AB",
            "primary_channel": "S_psi_psi^(N)/2=(S_AA+S_BB-2*S_AB)/2",
            "shell_definition": "equal integer squared norm |h|^2",
            "q_units": "nm^-1",
            "recorded_frames": (
                "direct_early_frames immediately preceding snapshots followed by "
                "direct_final_frames final snapshots (5+5 in the fixed study)"
            ),
            "convergence_frames": (
                "direct_early_frames immediately preceding the final window; final window alone "
                "defines the estimator"
            ),
            "common_q_bin_edges_inverse_nm": [
                0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.50
            ],
        },
        "runs": [asdict(spec) for spec in plan],
    }


def ensure_campaign(campaign_dir: Path, design: dict[str, Any], repo_root: Path) -> dict[str, Any]:
    manifest_path = campaign_dir / "manifest.json"
    design_hash = sha256_bytes(canonical_json(design))
    current_provenance = collect_provenance(repo_root)
    _validate_clean_git_provenance(current_provenance, "current")
    if manifest_path.exists():
        with manifest_path.open() as handle:
            existing = json.load(handle)
        try:
            validate_manifest_design_hash(existing)
        except ValueError as exc:
            raise RuntimeError(f"existing campaign manifest failed its design hash: {manifest_path}") from exc
        if existing.get("design_sha256") != design_hash:
            raise RuntimeError(
                f"campaign design differs from existing {manifest_path}; choose a new campaign id"
            )
        validate_creation_provenance(existing, current_provenance)
        return existing

    campaign_dir.mkdir(parents=True, exist_ok=True)
    (campaign_dir / "runs").mkdir(exist_ok=True)
    (campaign_dir / ".staging").mkdir(exist_ok=True)
    (campaign_dir / ".interrupted").mkdir(exist_ok=True)
    (campaign_dir / ".locks").mkdir(exist_ok=True)
    manifest = {
        **design,
        "design_sha256": design_hash,
        "created_utc": utc_now(),
        "creation_provenance": current_provenance,
        "creation_provenance_sha256": sha256_bytes(canonical_json(current_provenance)),
        "restart_semantics": (
            "completed run directories are checksum-verified and skipped; interrupted runs "
            "restart from the beginning when --recover-interrupted is supplied"
        ),
    }
    atomic_json(manifest_path, manifest)
    return manifest


def _archive_path(campaign_dir: Path, path: Path) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    destination = campaign_dir / ".interrupted" / f"{path.name}_{stamp}_{os.getpid()}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(destination))
    return destination


def expected_output_names(save_trajectory: bool) -> tuple[str, ...]:
    names = (
        "direct_structure_factor.npz",
        "meta.json",
        "planned_run.json",
        "snapshots.csv",
        "structure_factor.npz",
    )
    return (*names, "trajectory.npz") if save_trajectory else names


def completion_valid(
    run_dir: Path,
    expected_spec: RunSpec | dict[str, Any] | None = None,
    save_trajectory: bool | None = None,
    expected_git_commit: str | None = None,
) -> bool:
    completion_path = run_dir / "completion.json"
    if not completion_path.is_file():
        return False
    try:
        with completion_path.open() as handle:
            completion = json.load(handle)
        if completion.get("schema_version") != SCHEMA_VERSION:
            return False
        run_commit = _validate_clean_git_provenance(completion.get("provenance"), "run")
        if expected_git_commit is not None and run_commit != expected_git_commit:
            return False
        if not isinstance(completion.get("run_spec"), dict):
            return False
        if expected_spec is not None:
            expected = asdict(expected_spec) if isinstance(expected_spec, RunSpec) else expected_spec
            if completion.get("run_spec") != expected:
                return False
        recorded_save_trajectory = completion.get("save_trajectory")
        if not isinstance(recorded_save_trajectory, bool):
            return False
        if save_trajectory is not None and recorded_save_trajectory is not bool(save_trajectory):
            return False
        expected_names = expected_output_names(
            recorded_save_trajectory if save_trajectory is None else bool(save_trajectory)
        )
        if completion.get("required_outputs") != list(expected_names):
            return False
        files = completion["files"]
        if not isinstance(files, dict) or set(files) != set(expected_names):
            return False
        for relative, expected in files.items():
            if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
                return False
            if not isinstance(expected["bytes"], int) or expected["bytes"] <= 0:
                return False
            digest = expected["sha256"]
            if (
                not isinstance(digest, str)
                or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                return False
            path = run_dir / relative
            if not path.is_file() or sha256_file(path) != expected["sha256"]:
                return False
            if path.stat().st_size != expected["bytes"]:
                return False
    except (OSError, RuntimeError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return False
    return True


def acquire_lock(campaign_dir: Path, run_id: str, recover_interrupted: bool) -> Path:
    lock = campaign_dir / ".locks" / f"{run_id}.lock"
    if lock.exists():
        if not recover_interrupted:
            raise RuntimeError(
                f"lock exists for {run_id}: {lock}; verify no job is active, then use "
                "--recover-interrupted to restart it"
            )
        archived_lock = _archive_path(campaign_dir, lock)
        print(f"  archived stale lock -> {archived_lock}")
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError as exc:
        raise RuntimeError(f"another process acquired {lock}") from exc
    with os.fdopen(descriptor, "w") as handle:
        json.dump({"pid": os.getpid(), "hostname": socket.gethostname(), "created_utc": utc_now()}, handle)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    return lock


def required_output_paths(run_dir: Path, save_trajectory: bool) -> list[Path]:
    paths = [
        run_dir / "meta.json",
        run_dir / "snapshots.csv",
        run_dir / "structure_factor.npz",
        run_dir / "direct_structure_factor.npz",
        run_dir / "planned_run.json",
    ]
    if save_trajectory:
        paths.append(run_dir / "trajectory.npz")
    return paths


def run_spec(
    spec: RunSpec,
    args: argparse.Namespace,
    campaign_dir: Path,
    repo_root: Path,
) -> tuple[str, Path]:
    with (campaign_dir / "manifest.json").open() as handle:
        campaign_manifest = json.load(handle)
    validate_manifest_design_hash(campaign_manifest)
    campaign_commit = validate_creation_provenance(campaign_manifest)
    final_dir = campaign_dir / "runs" / spec.run_id
    if completion_valid(
        final_dir,
        spec,
        save_trajectory=bool(args.save_trajectory),
        expected_git_commit=campaign_commit,
    ):
        return "cached", final_dir
    if final_dir.exists() and not args.recover_interrupted:
        raise RuntimeError(
            f"incomplete or corrupt run exists at {final_dir}; use --recover-interrupted to archive and retry"
        )

    lock = acquire_lock(campaign_dir, spec.run_id, args.recover_interrupted)
    staging_dir = campaign_dir / ".staging" / spec.run_id
    try:
        if final_dir.exists():
            archived = _archive_path(campaign_dir, final_dir)
            print(f"  archived incomplete final directory -> {archived}")
        if staging_dir.exists():
            if not args.recover_interrupted:
                raise RuntimeError(
                    f"staging directory exists at {staging_dir}; use --recover-interrupted to archive and retry"
                )
            archived = _archive_path(campaign_dir, staging_dir)
            print(f"  archived interrupted staging directory -> {archived}")

        from .run import execute

        settings = simulation_settings(args)
        config = {
            **settings,
            "out": str(campaign_dir / ".staging"),
            "run_id": spec.run_id,
            "seed": spec.seed,
            "n_chains": spec.n_chains,
            "chain_length": spec.chain_length,
            "box_size": spec.box_size,
            "grid_size": spec.grid_size,
            "kappa": spec.kappa,
            "pi": spec.pi,
            "skip_if_cached": False,
        }
        staging_dir.mkdir(parents=True, exist_ok=False)
        started = utc_now()
        atomic_json(
            staging_dir / "planned_run.json",
            {"run_spec": asdict(spec), "simulation_settings": settings, "started_utc": started},
        )
        t0 = time.monotonic()
        returned = Path(execute(SimpleNamespace(**config))).resolve()
        if returned != staging_dir.resolve():
            raise RuntimeError(f"execute returned unexpected directory {returned}, expected {staging_dir}")

        outputs = required_output_paths(staging_dir, args.save_trajectory)
        missing = [str(path) for path in outputs if not path.is_file()]
        if missing:
            raise RuntimeError(f"run completed without required outputs: {missing}")
        file_manifest = {
            str(path.relative_to(staging_dir)): {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(outputs)
        }
        run_provenance = collect_provenance(repo_root)
        run_commit = _validate_clean_git_provenance(run_provenance, "run completion")
        if run_commit != campaign_commit:
            raise RuntimeError(
                f"run completed at commit {run_commit}, expected campaign commit {campaign_commit}"
            )
        completion = {
            "schema_version": SCHEMA_VERSION,
            "run_spec": asdict(spec),
            "save_trajectory": bool(args.save_trajectory),
            "required_outputs": list(expected_output_names(bool(args.save_trajectory))),
            "started_utc": started,
            "completed_utc": utc_now(),
            "elapsed_seconds": time.monotonic() - t0,
            "provenance": run_provenance,
            "files": file_manifest,
        }
        atomic_json(staging_dir / "completion.json", completion)
        final_dir.parent.mkdir(parents=True, exist_ok=True)
        os.replace(staging_dir, final_dir)
        return "completed", final_dir
    finally:
        try:
            lock.unlink()
        except FileNotFoundError:
            pass


def select_plan(plan: list[RunSpec], indices: list[int] | None) -> list[RunSpec]:
    if not indices:
        return plan
    invalid = [index for index in indices if index < 0 or index >= len(plan)]
    if invalid:
        raise ValueError(f"run indices out of range [0,{len(plan)-1}]: {invalid}")
    if len(set(indices)) != len(indices):
        raise ValueError("run indices must be unique")
    wanted = set(indices)
    return [spec for spec in plan if spec.index in wanted]


def print_plan(plan: list[RunSpec], design_hash: str) -> None:
    payload = {
        "design_sha256": design_hash,
        "n_runs": len(plan),
        "runs": [asdict(spec) for spec in plan],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Restartable fixed-density 144/288/576-chain finite-size study"
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="print the complete plan; write nothing")
    mode.add_argument(
        "--manifest-only", action="store_true", help="write/validate the campaign manifest; run nothing"
    )
    mode.add_argument("--run", action="store_true", help="execute selected production runs")
    mode.add_argument(
        "--analyze", action="store_true", help="aggregate completed direct spectra; run nothing"
    )
    parser.add_argument("--out", default="output/melt/fixed_density_size")
    parser.add_argument("--campaign-id", default="fixed_density_pi099")
    parser.add_argument("--sizes", type=int, nargs="+", default=list(DEFAULT_SIZES))
    parser.add_argument("--kappas", type=float, nargs="+", default=list(DEFAULT_KAPPAS))
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    parser.add_argument("--run-index", type=int, action="append", default=None)
    parser.add_argument("--chain-length", type=int, default=40)
    parser.add_argument("--reference-chains", type=int, default=144)
    parser.add_argument("--reference-box-size", type=float, default=22.0)
    parser.add_argument("--reference-grid-size", type=int, default=56)
    parser.add_argument("--pi", type=float, default=0.99)
    parser.add_argument("--f-A", dest="f_A", type=float, default=0.5)
    parser.add_argument("--block-length", type=int, default=4)
    parser.add_argument("--sequence-balance-max-attempts", type=int, default=100000)
    parser.add_argument("--bond-k", type=float, default=200.0)
    parser.add_argument("--bond-r0", type=float, default=1.0)
    parser.add_argument("--lj-eps-AA", dest="lj_eps_AA", type=float, default=1.0)
    parser.add_argument("--lj-eps-BB", dest="lj_eps_BB", type=float, default=1.0)
    parser.add_argument("--lj-eps-AB", dest="lj_eps_AB", type=float, default=0.1)
    parser.add_argument("--lj-eps-core", type=float, default=1.0)
    parser.add_argument("--lj-sigma", type=float, default=1.0)
    parser.add_argument("--lj-cutoff", type=float, default=2.5)
    parser.add_argument("--T-equilibrate", dest="T_equilibrate", type=float, default=5.0)
    parser.add_argument("--T-quench", dest="T_quench", type=float, default=0.7)
    parser.add_argument("--friction", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--equilibration", type=int, default=30000)
    parser.add_argument("--n-steps", type=int, default=250000)
    parser.add_argument("--snapshot-interval", type=int, default=2000)
    parser.add_argument("--platform", default=None)
    parser.add_argument("--direct-q-max", type=float, default=1.5)
    parser.add_argument("--direct-final-frames", type=int, default=5)
    parser.add_argument("--direct-early-frames", type=int, default=5)
    parser.add_argument("--direct-chunk-size", type=int, default=128)
    parser.add_argument("--save-trajectory", action="store_true")
    parser.add_argument("--recover-interrupted", action="store_true")
    parser.add_argument(
        "--allow-incomplete-analysis",
        action="store_true",
        help="with --analyze, aggregate available completed runs and record missing runs",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if float(args.f_A) != 0.5:
        raise ValueError("the fixed-density size study requires f_A=0.5")
    if float(args.lj_eps_AA) != float(args.lj_eps_BB):
        raise ValueError("the fixed-density size study requires lj_eps_AA=lj_eps_BB")
    if int(args.sequence_balance_max_attempts) < 1:
        raise ValueError("sequence_balance_max_attempts must be positive")
    if int(args.equilibration) < 0:
        raise ValueError("equilibration must be nonnegative")
    if int(args.n_steps) < 1 or int(args.snapshot_interval) < 1:
        raise ValueError("n_steps and snapshot_interval must be positive")
    if int(args.n_steps) % int(args.snapshot_interval) != 0:
        raise ValueError("n_steps must be exactly divisible by snapshot_interval")
    if args.direct_final_frames != 5 or args.direct_early_frames != 5:
        raise ValueError("the fixed-density study requires five early and five final direct frames")
    if args.n_steps // args.snapshot_interval < args.direct_final_frames + args.direct_early_frames:
        raise ValueError("production schedule has fewer snapshots than the two direct-frame windows")
    if args.direct_q_max <= 0.0 or args.direct_chunk_size < 1:
        raise ValueError("direct q maximum and chunk size must be positive")
    if float(args.direct_q_max) != 1.5:
        raise ValueError("the fixed-density study requires direct_q_max=1.5 nm^-1")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    validate_args(args)
    repo_root = Path(__file__).resolve().parents[1]
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
    design = campaign_design(args, plan)
    design_hash = sha256_bytes(canonical_json(design))
    selected = select_plan(plan, args.run_index)

    if args.dry_run:
        print_plan(plan, design_hash)
        return 0

    campaign_dir = Path(args.out).resolve() / args.campaign_id
    if args.analyze:
        manifest_path = campaign_dir / "manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"campaign manifest does not exist: {manifest_path}")
        from .fixed_density_analysis import aggregate_campaign

        outputs = aggregate_campaign(
            campaign_dir,
            require_complete=not args.allow_incomplete_analysis,
        )
        print(json.dumps({key: str(value) for key, value in outputs.items()}, indent=2, sort_keys=True))
        return 0

    manifest = ensure_campaign(campaign_dir, design, repo_root)
    print(
        f"campaign={args.campaign_id} design={manifest['design_sha256']} "
        f"runs={len(plan)} directory={campaign_dir}"
    )
    if args.manifest_only:
        print_plan(plan, design_hash)
        return 0
    print(f"selected_runs={len(selected)}")
    for position, spec in enumerate(selected, 1):
        print(
            f"[{position}/{len(selected)}] index={spec.index} {spec.run_id} "
            f"L={spec.box_size:.9f} G={spec.grid_size} rho={spec.bead_density:.12g}"
        )
        status, output = run_spec(spec, args, campaign_dir, repo_root)
        print(f"  {status}: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
