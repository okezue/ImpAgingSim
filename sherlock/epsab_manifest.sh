#!/bin/bash
# Step 1 of a campaign: write manifest.json once, before any array task starts.
# Light enough for a login node. Takes the same env vars as epsab_sweep.sbatch
# (CAMPAIGN_ID, EPS_ABS, SEEDS, SIZES, EXTRA_ARGS; see grid.sh).
#   bash sherlock/epsab_manifest.sh
#   CAMPAIGN_ID=epsab_stage2 EPS_ABS="0.5 0.475 0.45" SEEDS="1 2 3 4 5 6 7 8" SIZES="144 288" \
#     bash sherlock/epsab_manifest.sh
# Needs a clean git tree; the manifest pins the current commit.
set -euo pipefail

REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export REPO_DIR
source "$REPO_DIR/sherlock/grid.sh"
# No GPU here: record the platform the array tasks will use instead of probing.
export OPENMM_PLATFORM="${OPENMM_PLATFORM:-$CAMPAIGN_PLATFORM}"
source "$REPO_DIR/sherlock/env.sh"

if [ -n "$(git status --porcelain)" ]; then
  echo "epsab_manifest.sh: git tree is dirty; commit or remove these first:" >&2
  git status --short >&2
  exit 1
fi

grid_describe
set -x
python -m melt.epsab_scan --manifest-only \
  --out "$CAMPAIGN_ROOT" --campaign-id "$CAMPAIGN_ID" \
  --platform "$OPENMM_PLATFORM" \
  ${GRID_ARGS[@]+"${GRID_ARGS[@]}"} ${EXTRA_ARGS:-}
set +x

MANIFEST="$CAMPAIGN_ROOT/$CAMPAIGN_ID/manifest.json"
N_RUNS=$(python - "$MANIFEST" <<'EOF'
import json, sys
manifest = json.load(open(sys.argv[1]))
runs = manifest.get("runs")
print(len(runs) if isinstance(runs, list) else "?")
EOF
)
echo "manifest: $MANIFEST (runs=$N_RUNS)"
echo "next:     sbatch sherlock/epsab_sweep.sbatch        # same CAMPAIGN_ID/EPS_ABS/SEEDS/SIZES/EXTRA_ARGS"
