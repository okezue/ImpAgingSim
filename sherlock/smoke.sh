#!/bin/bash
# Smoke test on a compute node: unit tests, then a two-run campaign
# (eps_AB = 1.0 and 0.1, one seed, 4000 steps) under $CAMPAIGN_ROOT/smoke/smoke.
#   sh_dev -g 1                 # MIG GPU slice; plain sh_dev tests the CPU path
#   bash sherlock/smoke.sh
# Re-running is idempotent. After a code change the pinned commit no longer
# matches: rm -rf $SCRATCH/impagingsim/smoke, or set SMOKE_ID=smoke_<something>.
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export REPO_DIR
source "$REPO_DIR/sherlock/env.sh"

if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "smoke.sh: run this on a compute node (sh_dev -g 1), not a login node" >&2
  exit 2
fi

if [ -n "$(git status --porcelain)" ]; then
  echo "smoke.sh: git tree is dirty; the driver refuses to write a manifest. git status:" >&2
  git status --short >&2
  exit 1
fi

echo "=== pytest ==="
# -p no:cacheprovider keeps .pytest_cache out of the tree (it is not gitignored).
python -m pytest -p no:cacheprovider tests/ -q

SMOKE_ID="${SMOKE_ID:-smoke}"
SMOKE_OUT="$CAMPAIGN_ROOT/smoke"
SMOKE_ARGS=(
  --out "$SMOKE_OUT" --campaign-id "$SMOKE_ID"
  --eps-ABs 1.0 0.1 --seeds 1 --sizes 144
  --n-steps 4000 --snapshot-interval 1000 --mode-interval 200
  --platform "$OPENMM_PLATFORM"
)

echo "=== smoke campaign on $OPENMM_PLATFORM ==="
set -x
python -m melt.epsab_scan --manifest-only "${SMOKE_ARGS[@]}"
python -m melt.epsab_scan --run "${SMOKE_ARGS[@]}"
python -m melt.epsab_scan --analyze "${SMOKE_ARGS[@]}"
set +x

echo "=== outputs in $SMOKE_OUT/$SMOKE_ID ==="
find "$SMOKE_OUT/$SMOKE_ID" -type f -not -path '*/.*' | sort
du -sh "$SMOKE_OUT/$SMOKE_ID"
echo "smoke OK on platform $OPENMM_PLATFORM"
