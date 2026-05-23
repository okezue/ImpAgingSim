#!/bin/bash
# Corrected rerun campaign launched on AWS (g6e.2xlarge / L40S).
# Uses the fixed code path:
#   - reduced T* converted to Kelvin internally (T*=0.7 -> 84.19 K with eps_AA=1 kJ/mol)
#   - shared WCA repulsive core + tunable attractive tail (eps_AB only affects attraction)
#   - per-chain Markov sequence generation
#   - unified lambda=2*pi-1 generator across all f_A
#
# Each scan reproduces one decisive panel of the manuscript at corrected parameters.
# Launch on a single instance with:
#   GIT_REF=fix-rerun SCAN_KIND=rerun ./aws/launch.sh
# or run scans directly inside an already-launched instance.

set -euxo pipefail

PY=${PY:-/home/ubuntu/miniconda3/envs/imp/bin/python}
PLATFORM_FLAG=${PLATFORM_FLAG:---platform CUDA}
ROOT=${ROOT:-output/melt/scans_corrected}

mkdir -p ${ROOT}

# Fig 1: baseline kappa scan at corrected reduced units.
${PY} -m melt.kappa_scan \
  --scan_id fig1_baseline \
  --out ${ROOT} \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Fig 1 dense, soft, short-chain variants share the kappa axis with different geometry.
${PY} -m melt.kappa_scan \
  --scan_id fig1_dense \
  --out ${ROOT} \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 \
  --n_chains 144 --chain_length 40 --box_size 17.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

${PY} -m melt.kappa_scan \
  --scan_id fig1_soft \
  --out ${ROOT} \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AA 0.4 --lj_eps_BB 0.4 --lj_eps_AB 0.04 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

${PY} -m melt.kappa_scan \
  --scan_id fig1_short_chain \
  --out ${ROOT} \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 \
  --n_chains 480 --chain_length 12 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Fig 2: (pi, kappa) heatmap at f_A=0.5
${PY} -m melt.scan \
  --scan_id fig2_pi_kappa \
  --out ${ROOT} \
  --kind pi_kappa \
  --kappas 0.0 0.3 0.6 0.8 0.9 1.0 \
  --pis 0.5 0.7 0.85 0.95 0.98 0.99 0.995 0.999 \
  --seeds 1 2 3 4 5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Fig 3: (f_A, kappa) heatmap at pi=0.95
${PY} -m melt.scan \
  --scan_id fig3_fA_kappa \
  --out ${ROOT} \
  --kind fA_kappa \
  --kappas 0.0 0.2 0.4 0.6 0.8 1.0 \
  --f_As 0.2 0.3 0.4 0.5 0.6 0.7 0.8 \
  --pi 0.95 \
  --seeds 1 2 3 4 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Fig 4: epsAB scan with corrected WCA+tail force field. Lowering epsAB now reduces
# only the cross-attraction, not the cross-pair excluded volume.
${PY} -m melt.scan \
  --scan_id fig4_epsAB \
  --out ${ROOT} \
  --kind epsAB \
  --eps_ABs 0.0001 0.001 0.005 0.01 0.025 0.05 0.1 0.2 0.4 0.8 \
  --kappa 0.7 --pi 0.99 \
  --seeds 1 2 3 4 5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Fig 5 + dynamic observables: long aging with saved trajectories for F_s, MSD, chi_4.
# Save trajectories so melt.dynamics can compute F_s(k*,t;t_w), tau_alpha(t_w), MSD,
# alpha_2(t;t_w), chi_4 from saved positions. This addresses the reviewer's request
# for two-time dynamical observables before any glass-aging claim.
for kappa in 0.0 0.5 1.0; do
  for tw in 1000000 3000000 10000000 30000000 100000000; do
    for seed in 1 2 3; do
      ${PY} -m melt.run \
        --out ${ROOT}/fig5_aging \
        --run_id k${kappa}_pi0.99_tw${tw}_s${seed} \
        --sequence correlated --kappa ${kappa} --pi 0.99 --f_A 0.5 \
        --n_chains 144 --chain_length 40 --box_size 22.0 \
        --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
        --equilibration ${tw} --n_steps 250000 --snapshot_interval 2000 \
        --compute_density --grid_size 56 --save_trajectory \
        --seed ${seed} ${PLATFORM_FLAG}
    done
  done
done

aws s3 sync ${ROOT}/ "s3://${S3_BUCKET}/$(date +%Y-%m-%d)_corrected_$(hostname -s)/" --region "${REGION:-us-east-1}"
