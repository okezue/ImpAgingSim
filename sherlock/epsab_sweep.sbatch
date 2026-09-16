#!/bin/bash
#SBATCH --job-name=epsab
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --array=0-7
#SBATCH --output=slurm-%A_%a.out
# eps_AB sweep as a job array: task K runs plan indices i with i % S == K.
# epsab_manifest.sh must have been run first with the same env vars.
#   sbatch sherlock/epsab_sweep.sbatch
#   sbatch --array=0-15 --export=ALL,CAMPAIGN_ID=epsab_stage2,EPS_ABS="0.5 0.475 0.45",SEEDS="1 2 3 4 5 6 7 8",SIZES="144 288" \
#     sherlock/epsab_sweep.sbatch
#   sbatch --export=ALL,EXTRA_ARGS="--recover-interrupted" sherlock/epsab_sweep.sbatch   # after a timeout or crash
# Resubmitting is idempotent: completed runs are checksum-verified and skipped.
# Submit from the repo root; slurm-*.out logs land there and are gitignored.
set -euo pipefail

# $0 is a spool copy of this script under sbatch, so locate the kit via REPO_DIR.
REPO_DIR="${REPO_DIR:-$HOME/ImpAgingSim}"
export REPO_DIR
source "$REPO_DIR/sherlock/grid.sh"
source "$REPO_DIR/sherlock/env.sh"

MANIFEST="$CAMPAIGN_ROOT/$CAMPAIGN_ID/manifest.json"
if [ ! -f "$MANIFEST" ]; then
  echo "no manifest at $MANIFEST; run first: bash sherlock/epsab_manifest.sh (same CAMPAIGN_ID/EPS_ABS/SEEDS/SIZES)" >&2
  exit 1
fi
if [ "$OPENMM_PLATFORM" != "$CAMPAIGN_PLATFORM" ]; then
  echo "WARNING: running on $OPENMM_PLATFORM, not $CAMPAIGN_PLATFORM; check the CUDA install if this is a GPU node" >&2
fi

echo "job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} shard=${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT}" \
     "node=$(hostname -s) cpus=${SLURM_CPUS_PER_TASK:-?} start=$(date "+%Y-%m-%dT%H:%M:%S%z")"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader 2>/dev/null || echo "nvidia-smi unavailable"
grid_describe

set -x
python -m melt.epsab_scan --run \
  --out "$CAMPAIGN_ROOT" --campaign-id "$CAMPAIGN_ID" \
  --shard-index "$SLURM_ARRAY_TASK_ID" --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --platform "$OPENMM_PLATFORM" \
  ${GRID_ARGS[@]+"${GRID_ARGS[@]}"} ${EXTRA_ARGS:-}
set +x
echo "shard ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} done at $(date "+%Y-%m-%dT%H:%M:%S%z")"
