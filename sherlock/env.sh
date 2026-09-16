#!/bin/bash
# Sourced by every Sherlock job and interactive shell: python module, venv,
# OpenMM platform, campaign paths, cd into the repo.
#   source sherlock/env.sh
#
#   REPO_DIR         default $HOME/ImpAgingSim
#   VENV             default $HOME/venvs/imp (built by setup_env.sh)
#   PY_MODULE        default: module recorded by setup_env.sh, else python/3.12.1
#   CAMPAIGN_ROOT    default $SCRATCH/impagingsim; must stay outside the git tree
#   OPENMM_PLATFORM  preset to skip the probe; else CUDA -> OpenCL -> CPU
case $- in *i*) ;; *) set -euo pipefail ;; esac

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  echo "env.sh is meant to be sourced: source sherlock/env.sh" >&2
  exit 2
fi

_env_fail() {
  echo "env.sh: $*" >&2
  return 1
}

# Lmod's module function and venv activate are not reliably `set -u` clean.
_env_lenient() {
  local opts=$-
  set +u
  "$@"
  local rc=$?
  case $opts in *u*) set -u ;; esac
  return $rc
}

export REPO_DIR="${REPO_DIR:-$HOME/ImpAgingSim}"
export VENV="${VENV:-$HOME/venvs/imp}"
[ -n "${SCRATCH:-}" ] || _env_fail "\$SCRATCH is not set; this file is for Sherlock" || return 1
export CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-$SCRATCH/impagingsim}"
case "$CAMPAIGN_ROOT" in
  "$REPO_DIR"/*) _env_fail "CAMPAIGN_ROOT=$CAMPAIGN_ROOT is inside the repo; outputs would dirty the git tree" || return 1 ;;
esac

if ! type module >/dev/null 2>&1; then
  for _env_init in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do
    if [ -r "$_env_init" ]; then
      _env_lenient source "$_env_init"
      break
    fi
  done
  unset _env_init
fi
type module >/dev/null 2>&1 || _env_fail "'module' not found; start from a login shell (bash -l)" || return 1

if [ -z "${PY_MODULE:-}" ]; then
  if [ -r "$VENV/sherlock_python_module" ]; then
    PY_MODULE=$(<"$VENV/sherlock_python_module")
  else
    PY_MODULE=python/3.12.1
  fi
fi
export PY_MODULE
_env_lenient module load "$PY_MODULE" || _env_fail "cannot load $PY_MODULE (see: ml spider python)" || return 1

[ -f "$VENV/bin/activate" ] || _env_fail "no venv at $VENV; run 'bash sherlock/setup_env.sh' on an sh_dev node" || return 1
_env_lenient source "$VENV/bin/activate"

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg
# Do not oversubscribe shared nodes: numpy/FFT and the OpenMM CPU platform.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export OPENMM_CPU_THREADS="${OPENMM_CPU_THREADS:-$OMP_NUM_THREADS}"
if [ -n "${L_SCRATCH:-}" ] && [ -d "$L_SCRATCH" ]; then
  export TMPDIR="$L_SCRATCH"
fi

# Tiny Context on the named platform; prints "<name> OK (<device>)" on success.
_env_probe_platform() {
  local runner=(python -)
  if command -v timeout >/dev/null 2>&1; then
    runner=(timeout "${OPENMM_PROBE_TIMEOUT:-300}" python -)
  fi
  "${runner[@]}" "$1" <<'EOF' 2>/dev/null
import sys
import openmm as mm

name = sys.argv[1]
system = mm.System()
system.addParticle(1.0)
integrator = mm.LangevinMiddleIntegrator(1.0, 1.0, 0.005)
context = mm.Context(system, integrator, mm.Platform.getPlatformByName(name))
device = "host"
if name in ("CUDA", "OpenCL"):
    device = context.getPlatform().getPropertyValue(context, "DeviceName")
print(f"{name} OK ({device})")
EOF
}

if [ -z "${OPENMM_PLATFORM:-}" ]; then
  for _env_candidate in CUDA OpenCL CPU; do
    if _env_probe_platform "$_env_candidate"; then
      OPENMM_PLATFORM=$_env_candidate
      break
    fi
  done
  unset _env_candidate
  if [ -z "${OPENMM_PLATFORM:-}" ]; then
    echo "env.sh: WARNING no CUDA/OpenCL/CPU platform works, falling back to Reference" >&2
    OPENMM_PLATFORM=Reference
  fi
fi
export OPENMM_PLATFORM

[ -d "$REPO_DIR/melt" ] || _env_fail "REPO_DIR=$REPO_DIR does not contain melt/" || return 1
cd "$REPO_DIR"
mkdir -p "$CAMPAIGN_ROOT"

echo "env: $PY_MODULE $(python --version 2>&1) venv=$VENV platform=$OPENMM_PLATFORM" \
  "repo=$REPO_DIR@$(git rev-parse --short HEAD 2>/dev/null || echo '?') campaign_root=$CAMPAIGN_ROOT"
