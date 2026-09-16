#!/bin/bash
# One-time, idempotent environment build. Run on a compute node:
#   sh_dev                      # or: sh_dev -g 1  to also exercise the CUDA platform
#   bash sherlock/setup_env.sh
# Loads a python module (3.10 <= v <= 3.12: numpy<2 needs <=3.12, OpenMM PyPI
# wheels need >=3.10), creates ${VENV:-$HOME/venvs/imp}, installs
# requirements.txt + openmm[cuda12] + pytest, then prints the OpenMM self test.
#   PY_MODULE=python/3.11.x     force a module; otherwise python/3.12.1 then ml spider
#   FORCE_LOGIN_NODE=1          allow running outside a SLURM job
set -euo pipefail

if [ -z "${SLURM_JOB_ID:-}" ] && [ "${FORCE_LOGIN_NODE:-0}" != "1" ]; then
  cat >&2 <<'EOF'
setup_env.sh: not inside a SLURM job. Sherlock asks that environments be built
on a compute node, not a login node:
    sh_dev
    bash sherlock/setup_env.sh
Set FORCE_LOGIN_NODE=1 to override.
EOF
  exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR="${REPO_DIR:-$(dirname "$SCRIPT_DIR")}"
VENV="${VENV:-$HOME/venvs/imp}"
# Keep the wheel cache (CUDA runtime wheels are large) off the 15 GB $HOME.
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${SCRATCH:-$HOME}/.cache/pip}"

[ -f "$REPO_DIR/requirements.txt" ] || { echo "no requirements.txt in REPO_DIR=$REPO_DIR" >&2; exit 1; }
if [ "$REPO_DIR" != "$HOME/ImpAgingSim" ]; then
  echo "NOTE: jobs default to REPO_DIR=\$HOME/ImpAgingSim; this clone is $REPO_DIR." \
       "Add 'export REPO_DIR=$REPO_DIR' to ~/.bashrc or pass it with --export."
fi

# Lmod's module function and venv activate are not reliably `set -u` clean.
lenient() {
  local opts=$-
  set +u
  "$@"
  local rc=$?
  case $opts in *u*) set -u ;; esac
  return $rc
}

if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do
    if [ -r "$init" ]; then lenient source "$init"; break; fi
  done
fi
type module >/dev/null 2>&1 || { echo "'module' not found; start from a login shell (bash -l)" >&2; exit 1; }

# python/3.x.y with 10 <= x <= 12
version_ok() {
  local minor
  minor=$(sed -nE 's#^python/3\.([0-9]+)(\..*)?$#\1#p' <<<"$1")
  [ -n "$minor" ] && [ "$minor" -ge 10 ] && [ "$minor" -le 12 ]
}

candidates=()
if [ -n "${PY_MODULE:-}" ]; then candidates+=("$PY_MODULE"); fi
# An existing venv must keep the interpreter it was built with.
if [ -r "$VENV/sherlock_python_module" ]; then candidates+=("$(<"$VENV/sherlock_python_module")"); fi
candidates+=(python/3.12.1)
while IFS= read -r m; do
  if [ -n "$m" ]; then candidates+=("$m"); fi
done < <(lenient module -t spider python 2>&1 | grep -oE 'python/3\.[0-9]+(\.[0-9]+)?' | sort -urV || true)

PY_MODULE_LOADED=""
for m in "${candidates[@]}"; do
  if ! version_ok "$m"; then
    echo "skipping $m (need python 3.10 - 3.12)"
    continue
  fi
  if lenient module load "$m" 2>/dev/null; then
    PY_MODULE_LOADED=$m
    break
  fi
  echo "module load $m failed, trying next"
done
[ -n "$PY_MODULE_LOADED" ] || { echo "no loadable python module in 3.10 - 3.12; check: ml spider python" >&2; exit 1; }
echo "=== using $PY_MODULE_LOADED: $(command -v python3) ($(python3 --version 2>&1)) ==="

if [ -f "$VENV/bin/activate" ]; then
  echo "=== venv exists at $VENV, reusing ==="
else
  echo "=== creating venv at $VENV ==="
  mkdir -p "$(dirname "$VENV")"
  python3 -m venv "$VENV"
fi
lenient source "$VENV/bin/activate"

echo "=== installing packages ==="
python -m pip install --upgrade pip
python -m pip install -r "$REPO_DIR/requirements.txt" "openmm[cuda12]" pytest
echo "$PY_MODULE_LOADED" > "$VENV/sherlock_python_module"

echo "=== nvidia-smi ==="
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv 2>/dev/null || echo "(no GPU on this node)"
echo "=== python -m openmm.testInstallation ==="
python -m openmm.testInstallation || echo "WARNING: testInstallation reported errors (CUDA/OpenCL failures are expected on a node without a GPU)"
echo "=== platforms ==="
python -c "import openmm as mm; print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])"
echo "=== venv size ==="
du -sh "$VENV"
echo "next: sh_dev -g 1 ; bash sherlock/smoke.sh"
