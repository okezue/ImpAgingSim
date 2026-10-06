#!/bin/bash
# Static parameter sweeps supporting the manuscript and Supplementary Information.
# Uses the fixed code path:
#   - reduced T* converted to Kelvin internally (T*=0.7 -> 84.19 K with eps_AA=1 kJ/mol)
#   - shared WCA repulsive core + tunable attractive tail (eps_AB only affects attraction)
#   - per-chain Markov sequence generation
#   - unified lambda=2*pi-1 generator across all f_A
#
# Historical scan IDs are retained so the source-data scripts can locate results.
# See aws/RESULTS_MANIFEST.md for their current manuscript figure mapping.
# Run locally in the simulation environment or use SCAN_KIND=rerun on EC2.

set -euxo pipefail

PY=${PY:-python}
PLATFORM_FLAG=${PLATFORM_FLAG:---platform CUDA}
ROOT=${ROOT:-output/melt/scans_corrected}

mkdir -p "${ROOT}"

# Main Fig. 2: baseline kappa scan at corrected reduced units.
${PY} -m melt.kappa_scan \
  --scan_id fig1_baseline \
  --out "${ROOT}" \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --pi 0.90 --f_A 0.5 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Main Fig. 2 and Supplementary Fig. S5: dense, soft and short-chain controls.
${PY} -m melt.kappa_scan \
  --scan_id fig1_dense \
  --out "${ROOT}" \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 \
  --n_chains 144 --chain_length 40 --box_size 17.0 \
  --pi 0.90 --f_A 0.5 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# The current driver references the unit-depth core; T*=0.28 preserves 33.676 K.
${PY} -m melt.kappa_scan \
  --scan_id fig1_soft \
  --out "${ROOT}" \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --pi 0.90 --f_A 0.5 \
  --T_equilibrate 2.0 --T_quench 0.28 --lj_eps_AA 0.4 --lj_eps_BB 0.4 --lj_eps_AB 0.04 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

${PY} -m melt.kappa_scan \
  --scan_id fig1_short_chain \
  --out "${ROOT}" \
  --kappas 0.0 0.15 0.3 0.45 0.6 0.75 0.9 1.0 \
  --seeds 1 2 3 4 \
  --n_chains 480 --chain_length 12 --box_size 22.0 \
  --pi 0.90 --f_A 0.5 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Main Figs. 3-4: (pi, kappa) grid at f_A=0.5
${PY} -m melt.scan \
  --scan_id fig2_pi_kappa \
  --out "${ROOT}" \
  --kind pi_kappa \
  --kappas 0.0 0.3 0.6 0.8 0.9 1.0 \
  --pis 0.5 0.7 0.85 0.95 0.98 0.99 0.995 0.999 \
  --seeds 1 2 3 4 5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Supplementary Fig. S6: off-stoichiometric mode-mixing diagnostic at pi=0.95
${PY} -m melt.scan \
  --scan_id fig3_fA_kappa \
  --out "${ROOT}" \
  --kind fA_kappa \
  --kappas 0.0 0.2 0.4 0.6 0.8 1.0 \
  --f_As 0.2 0.3 0.4 0.5 0.6 0.7 0.8 \
  --pi 0.95 \
  --seeds 1 2 3 4 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

# Main Fig. 5 and Supplementary Fig. S7: sequence-class / epsAB scan. Lowering epsAB reduces
# only the cross-attraction, not the cross-pair excluded volume.
${PY} -m melt.scan \
  --scan_id fig4_epsAB \
  --out "${ROOT}" \
  --kind epsAB \
  --eps_ABs 0.0001 0.001 0.005 0.01 0.025 0.05 0.1 0.2 0.4 0.8 \
  --sequences alternating random block correlated \
  --kappa 0.7 --pi 0.99 \
  --seeds 1 2 3 4 5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 \
  --n_steps 250000 --equilibration 30000 --snapshot_interval 2000 \
  --grid_size 56 ${PLATFORM_FLAG}

if [ -n "${S3_BUCKET:-}" ]; then
  aws s3 sync "${ROOT}/" "s3://${S3_BUCKET}/$(date +%Y-%m-%d)_corrected_$(hostname -s)/" --region "${REGION:-us-east-1}"
fi
