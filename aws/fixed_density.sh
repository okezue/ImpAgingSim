#!/bin/bash
# Fixed-density finite-size study: M = 144/288/576 chains at constant bead
# density, kappa = 0 and 1, five seeds each. 30 production runs.
#
# Box scales as L(M) = 22 (M/144)^(1/3) and the density grid scales with it, so
# bead density and real-space mesh spacing are held fixed across system sizes.
#
# Restartable: runs are staged, hashed and atomically promoted, so re-invoking
# this script verifies hashes and skips completed runs.

set -euxo pipefail

PY=${PY:-/home/ubuntu/miniconda3/envs/imp/bin/python}
PLATFORM_FLAG=${PLATFORM_FLAG:---platform CUDA}
CAMPAIGN=${CAMPAIGN:-fixed_density_pi099_v2}

${PY} -m melt.fixed_density_size_scan --manifest-only --campaign-id ${CAMPAIGN} ${PLATFORM_FLAG}
${PY} -m melt.fixed_density_size_scan --run --campaign-id ${CAMPAIGN} ${PLATFORM_FLAG}
${PY} -m melt.fixed_density_size_scan --analyze --campaign-id ${CAMPAIGN} ${PLATFORM_FLAG}

aws s3 sync output/melt/fixed_density_size/ \
  "s3://${S3_BUCKET}/$(date +%Y-%m-%d)_fixed_density_$(hostname -s)/" \
  --region "${REGION:-us-east-1}"
