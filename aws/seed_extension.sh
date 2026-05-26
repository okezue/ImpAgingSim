#!/bin/bash
# Seed extension: add seeds 4-10 at t_w=10^8 for each kappa in {0.0, 0.5, 1.0}
# to tighten the chi_4^peak amplitude estimate (currently from 3 seeds).
#
# 3 kappa x 7 seeds = 21 runs. Each is t_w=10^8 BD equilibration + 250k production,
# roughly 3.5 hr per run on g6e.2xlarge L40S. ~75 GPU-hr total.
#
# Idempotent: uses melt.run skip_if_cached, so re-launching after a crash skips
# anything already in output/melt/scans_corrected/fig5_aging/.

set -euxo pipefail

PY=${PY:-/home/ubuntu/miniconda3/envs/imp/bin/python}
PLATFORM_FLAG=${PLATFORM_FLAG:---platform CUDA}
ROOT=${ROOT:-output/melt/scans_corrected}

mkdir -p ${ROOT}/fig5_aging

for kappa in 0.0 0.5 1.0; do
  for seed in 4 5 6 7 8 9 10; do
    ${PY} -m melt.run \
      --out ${ROOT}/fig5_aging \
      --run_id k${kappa}_pi0.99_tw100000000_s${seed} \
      --sequence correlated --kappa ${kappa} --pi 0.99 --f_A 0.5 \
      --n_chains 144 --chain_length 40 --box_size 22.0 \
      --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
      --equilibration 100000000 --n_steps 250000 --snapshot_interval 2000 \
      --compute_density --grid_size 56 --save_trajectory \
      --seed ${seed} ${PLATFORM_FLAG}
  done
done

aws s3 sync ${ROOT}/fig5_aging/ "s3://${S3_BUCKET}/$(date +%Y-%m-%d)_seed_extension_$(hostname -s)/" --region "${REGION:-us-east-1}"
