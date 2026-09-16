#!/bin/bash
#SBATCH --job-name=epsab_analyze
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --output=slurm-%j.out
# Aggregate a finished (or partially finished) eps_AB campaign. CPU only.
#   sbatch sherlock/epsab_analyze.sbatch
#   sbatch --dependency=afterok:<array job id> sherlock/epsab_analyze.sbatch
#   sbatch --export=ALL,ALLOW_INCOMPLETE=1 sherlock/epsab_analyze.sbatch      # aggregate what is complete
# Same CAMPAIGN_ID/EPS_ABS/SEEDS/SIZES/EXTRA_ARGS as the manifest (see grid.sh).
set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/ImpAgingSim}"
export REPO_DIR
source "$REPO_DIR/sherlock/grid.sh"
# Analysis runs no dynamics; the platform flag is accepted but unused here.
export OPENMM_PLATFORM="${OPENMM_PLATFORM:-$CAMPAIGN_PLATFORM}"
source "$REPO_DIR/sherlock/env.sh"

CAMPAIGN_DIR="$CAMPAIGN_ROOT/$CAMPAIGN_ID"
if [ ! -f "$CAMPAIGN_DIR/manifest.json" ]; then
  echo "no manifest at $CAMPAIGN_DIR/manifest.json; nothing to analyze" >&2
  exit 1
fi

ANALYZE_ARGS=()
if [ "${ALLOW_INCOMPLETE:-0}" = "1" ]; then
  ANALYZE_ARGS+=(--allow-incomplete-analysis)
fi

echo "job=${SLURM_JOB_ID} node=$(hostname -s) start=$(date "+%Y-%m-%dT%H:%M:%S%z")"
grid_describe
echo "completed runs: $(find "$CAMPAIGN_DIR/runs" -mindepth 2 -maxdepth 2 -name completion.json 2>/dev/null | wc -l | tr -d " ")"

set -x
python -m melt.epsab_scan --analyze \
  --out "$CAMPAIGN_ROOT" --campaign-id "$CAMPAIGN_ID" \
  --platform "$OPENMM_PLATFORM" \
  ${GRID_ARGS[@]+"${GRID_ARGS[@]}"} ${ANALYZE_ARGS[@]+"${ANALYZE_ARGS[@]}"} ${EXTRA_ARGS:-}
set +x

echo "=== $CAMPAIGN_DIR/analysis ==="
ls -la "$CAMPAIGN_DIR/analysis" 2>/dev/null || echo "(no analysis/ directory; see the driver output above for where it wrote)"
echo "analysis done at $(date "+%Y-%m-%dT%H:%M:%S%z")"
