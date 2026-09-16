"""Generic restartable-campaign machinery shared by the scan drivers.

A campaign directory holds ``manifest.json`` (the hashed design plus the git commit it was
created at), ``runs/<run_id>/`` for completed runs, and ``.staging``, ``.locks``,
``.interrupted`` for in-flight work.  Every run is executed in staging, its required
outputs are hashed into ``completion.json``, and the directory is atomically promoted.
Re-invoking a campaign verifies hashes and skips finished runs, so shards of a scheduler
array or several processes on one filesystem can share a campaign safely.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import time
from types import SimpleNamespace
from typing import Any, Callable

from .fixed_density_size_scan import (
    _archive_path,
    _validate_clean_git_provenance,
    acquire_lock,
    atomic_json,
    canonical_json,
    collect_provenance,
    sha256_bytes,
    sha256_file,
    utc_now,
    validate_creation_provenance,
)

SCHEMA_VERSION = 1


def design_hash(design: dict[str, Any], keys: tuple[str, ...]) -> str:
    missing = [key for key in keys if key not in design]
    if missing:
        raise ValueError(f"campaign design lacks fields {missing}")
    return sha256_bytes(canonical_json({key: design[key] for key in keys}))


def validate_manifest(manifest: dict[str, Any], keys: tuple[str, ...]) -> str:
    if not isinstance(manifest, dict):
        raise ValueError("campaign manifest root must be an object")
    actual = design_hash(manifest, keys)
    if manifest.get("design_sha256") != actual:
        raise ValueError("campaign manifest failed its design hash")
    return actual


def ensure_campaign(
    campaign_dir: Path,
    design: dict[str, Any],
    keys: tuple[str, ...],
    repo_root: Path,
    require_clean_git: bool = True,
) -> dict[str, Any]:
    manifest_path = campaign_dir / "manifest.json"
    digest = design_hash(design, keys)
    current = collect_provenance(repo_root)
    if require_clean_git:
        _validate_clean_git_provenance(current, "current")
    if manifest_path.exists():
        with manifest_path.open() as handle:
            existing = json.load(handle)
        try:
            validate_manifest(existing, keys)
        except ValueError as exc:
            raise RuntimeError(f"existing campaign manifest failed its design hash: {manifest_path}") from exc
        if existing.get("design_sha256") != digest:
            raise RuntimeError(
                f"campaign design differs from existing {manifest_path}; choose a new campaign id"
            )
        if require_clean_git:
            validate_creation_provenance(existing, current)
        return existing
    for name in ("runs", ".staging", ".interrupted", ".locks"):
        (campaign_dir / name).mkdir(parents=True, exist_ok=True)
    manifest = {
        **design,
        "design_sha256": digest,
        "created_utc": utc_now(),
        "creation_provenance": current,
        "creation_provenance_sha256": sha256_bytes(canonical_json(current)),
        "require_clean_git": bool(require_clean_git),
        "restart_semantics": (
            "completed run directories are checksum-verified and skipped; interrupted runs "
            "restart from the beginning when --recover-interrupted is supplied"
        ),
    }
    atomic_json(manifest_path, manifest)
    return manifest


def campaign_commit(manifest: dict[str, Any]) -> str | None:
    """Commit pinned at creation, or None for campaigns created without the clean-tree guard."""
    if manifest.get("require_clean_git", True):
        return validate_creation_provenance(manifest)
    return None


def completion_valid(
    run_dir: Path,
    expected_spec: dict[str, Any] | None,
    required_outputs: tuple[str, ...],
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
        if expected_git_commit is not None:
            run_commit = _validate_clean_git_provenance(completion.get("provenance"), "run")
            if run_commit != expected_git_commit:
                return False
        if expected_spec is not None and completion.get("run_spec") != expected_spec:
            return False
        if completion.get("required_outputs") != list(required_outputs):
            return False
        files = completion["files"]
        if not isinstance(files, dict) or set(files) != set(required_outputs):
            return False
        for relative, expected in files.items():
            path = run_dir / relative
            if not path.is_file() or path.stat().st_size != expected["bytes"]:
                return False
            if sha256_file(path) != expected["sha256"]:
                return False
    except (OSError, RuntimeError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return False
    return True


def run_staged(
    run_id: str,
    spec: dict[str, Any],
    config: dict[str, Any],
    required_outputs: tuple[str, ...],
    campaign_dir: Path,
    repo_root: Path,
    expected_commit: str | None,
    recover_interrupted: bool,
    execute: Callable[[Any], str] | None = None,
) -> tuple[str, Path]:
    """Run one spec in staging and promote it atomically; return ("cached"|"completed", dir)."""
    final_dir = campaign_dir / "runs" / run_id
    if completion_valid(final_dir, spec, required_outputs, expected_commit):
        return "cached", final_dir
    if final_dir.exists() and not recover_interrupted:
        raise RuntimeError(
            f"incomplete or corrupt run exists at {final_dir}; use --recover-interrupted to archive and retry"
        )
    if execute is None:
        from .run import execute as execute_run

        execute = execute_run
    lock = acquire_lock(campaign_dir, run_id, recover_interrupted)
    staging_dir = campaign_dir / ".staging" / run_id
    try:
        if final_dir.exists():
            print(f"  archived incomplete final directory -> {_archive_path(campaign_dir, final_dir)}")
        if staging_dir.exists():
            if not recover_interrupted:
                raise RuntimeError(
                    f"staging directory exists at {staging_dir}; use --recover-interrupted to archive and retry"
                )
            print(f"  archived interrupted staging directory -> {_archive_path(campaign_dir, staging_dir)}")
        staging_dir.mkdir(parents=True, exist_ok=False)
        started = utc_now()
        atomic_json(
            staging_dir / "planned_run.json",
            {"run_spec": spec, "config": _jsonable(config), "started_utc": started},
        )
        t0 = time.monotonic()
        returned = Path(execute(SimpleNamespace(**config, out=str(campaign_dir / ".staging"), run_id=run_id))).resolve()
        if returned != staging_dir.resolve():
            raise RuntimeError(f"execute returned unexpected directory {returned}, expected {staging_dir}")
        missing = [name for name in required_outputs if not (staging_dir / name).is_file()]
        if missing:
            raise RuntimeError(f"run completed without required outputs: {missing}")
        provenance = collect_provenance(repo_root)
        if expected_commit is not None:
            run_commit = _validate_clean_git_provenance(provenance, "run completion")
            if run_commit != expected_commit:
                raise RuntimeError(f"run completed at commit {run_commit}, expected {expected_commit}")
        completion = {
            "schema_version": SCHEMA_VERSION,
            "run_spec": spec,
            "required_outputs": list(required_outputs),
            "started_utc": started,
            "completed_utc": utc_now(),
            "elapsed_seconds": time.monotonic() - t0,
            "provenance": provenance,
            "files": {
                name: {"bytes": (staging_dir / name).stat().st_size, "sha256": sha256_file(staging_dir / name)}
                for name in required_outputs
            },
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


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def select_shard(indices: list[int], shard_index: int | None, shard_count: int | None) -> list[int]:
    if shard_index is None and shard_count is None:
        return indices
    if shard_index is None or shard_count is None:
        raise ValueError("--shard-index and --shard-count must be given together")
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError(f"shard index {shard_index} must lie in [0,{shard_count})")
    return [index for index in indices if index % shard_count == shard_index]
