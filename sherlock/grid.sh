#!/bin/bash
# Shared by epsab_manifest.sh, epsab_sweep.sbatch and epsab_analyze.sbatch.
# Turns the campaign env vars into melt.epsab_scan flags so every phase sees
# the same design (designs are hashed; the driver rejects a mismatch).
#
#   CAMPAIGN_ID         default epsab_kappa05
#   CAMPAIGN_PLATFORM   --platform passed where no GPU can be probed (manifest,
#                       analyze); default CUDA. Not part of the design hash.
#   EPS_ABS SEEDS SIZES space separated lists; unset -> driver defaults
#   EXTRA_ARGS          appended verbatim to every phase, e.g. "--n-steps 500000"
#
# Source it, do not execute it. Sets CAMPAIGN_ID, CAMPAIGN_PLATFORM, GRID_ARGS.
case $- in *i*) ;; *) set -euo pipefail ;; esac

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  echo "grid.sh is meant to be sourced: source sherlock/grid.sh" >&2
  exit 2
fi

export CAMPAIGN_ID="${CAMPAIGN_ID:-epsab_kappa05}"
export CAMPAIGN_PLATFORM="${CAMPAIGN_PLATFORM:-CUDA}"

GRID_ARGS=()
if [ -n "${EPS_ABS:-}" ]; then
  read -r -a _grid_eps <<<"${EPS_ABS}"
  GRID_ARGS+=(--eps-ABs "${_grid_eps[@]}")
fi
if [ -n "${SEEDS:-}" ]; then
  read -r -a _grid_seeds <<<"${SEEDS}"
  GRID_ARGS+=(--seeds "${_grid_seeds[@]}")
fi
if [ -n "${SIZES:-}" ]; then
  read -r -a _grid_sizes <<<"${SIZES}"
  GRID_ARGS+=(--sizes "${_grid_sizes[@]}")
fi
unset _grid_eps _grid_seeds _grid_sizes

grid_describe() {
  echo "campaign_id=${CAMPAIGN_ID} manifest_platform=${CAMPAIGN_PLATFORM}"
  echo "eps_ABs=${EPS_ABS:-<driver default: 37 points 1.0 -> 0.1>}"
  echo "seeds=${SEEDS:-<driver default: 1 2 3 4>}"
  echo "sizes=${SIZES:-<driver default: 144>}"
  echo "extra_args=${EXTRA_ARGS:-<none>}"
}
