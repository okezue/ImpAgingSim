#!/bin/bash
# Run on the LAPTOP: pull a campaign from Sherlock $SCRATCH through the data
# transfer nodes into output/sherlock/<campaign>/.
#   SUNETID=<sunetid> bash sherlock/sync_back.sh            # small files only (default)
#   SUNETID=<sunetid> FULL=1 bash sherlock/sync_back.sh     # everything, incl. mode_amplitudes.npz (~3-4 GB for stage 1)
#   CAMPAIGN=epsab_stage2 ...                               # default epsab_kappa05
# Small files = manifest.json, analysis/, runs/*/{completion.json,meta.json,snapshots.csv}.
# $SCRATCH on Sherlock is /scratch/users/<sunetid>; set REMOTE_ROOT if CAMPAIGN_ROOT was changed.
# Each invocation is one Duo prompt unless ~/.ssh/config multiplexes the DTN (see README).
set -euo pipefail

SUNETID="${SUNETID:?must export SUNETID=<your sunetid>}"
CAMPAIGN="${CAMPAIGN:-epsab_kappa05}"
DTN="${DTN:-dtn.sherlock.stanford.edu}"
REMOTE_ROOT="${REMOTE_ROOT:-/scratch/users/${SUNETID}/impagingsim}"
LOCAL_DIR="${LOCAL_DIR:-output/sherlock/${CAMPAIGN}}"
SRC="${SUNETID}@${DTN}:${REMOTE_ROOT}/${CAMPAIGN}/"

FILTERS=(--exclude='.staging/' --exclude='.locks/' --exclude='.interrupted/')
if [ "${FULL:-0}" = "1" ]; then
  MODE="full"
else
  MODE="small files only (FULL=1 for everything)"
  FILTERS+=(
    --include='manifest.json'
    --include='analysis/***'
    --include='runs/'
    --include='runs/*/'
    --include='runs/*/completion.json'
    --include='runs/*/meta.json'
    --include='runs/*/snapshots.csv'
    --exclude='*'
    --prune-empty-dirs
  )
fi

mkdir -p "${LOCAL_DIR}"
echo "Syncing ${SRC} -> ${LOCAL_DIR}/ [${MODE}] ..."
rsync -avP "${FILTERS[@]}" "${SRC}" "${LOCAL_DIR}/"
echo "Done."
du -sh "${LOCAL_DIR}"
echo "completed runs pulled: $(find "${LOCAL_DIR}/runs" -name completion.json 2>/dev/null | wc -l | tr -d ' ')"
