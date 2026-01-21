#!/bin/bash
# Robustness simulation runner for IMP aging study
# This script runs all the simulations specified in the dev spec

set -e

OUTDIR="output/robustness"
mkdir -p "$OUTDIR"

echo "=============================================="
echo "IMP Aging Robustness Simulations"
echo "=============================================="
echo ""

# Common parameters
SEED=12345
S=20        # n_disorder
M=8         # n_traj
PRE=2000    # pre_sweeps
BETA0=0.05
BETA=1.0

# ==============================================================================
# A2: Finite-size robustness (N = 20, 30, 40)
# ==============================================================================
echo ">>> A2: Finite-size robustness (N=20,30,40)"
echo ""

for N in 20 30 40; do
    for EPS in 3 6; do
        for ENS in iid correlated; do
            OUTFILE="${OUTDIR}/results_N${N}_eps${EPS}_${ENS}.csv"
            echo "Running: N=$N, eps=$EPS, ens=$ENS -> $OUTFILE"
            python3 imp_aging.py \
                --N $N \
                --ensemble $ENS \
                --epsilon $EPS \
                --n_disorder $S \
                --n_traj $M \
                --pre_sweeps $PRE \
                --beta0 $BETA0 \
                --beta $BETA \
                --tw 0 100 300 1000 3000 10000 \
                --t_max 10000 \
                --n_lags 28 \
                --kappa 0.7 \
                --pi 0.9 \
                --seed $SEED \
                --out "$OUTFILE"
        done
    done
done

echo ""
echo "A2 complete!"
echo ""

# ==============================================================================
# A3: Correlated-disorder parameter scan (kappa, pi)
# ==============================================================================
echo ">>> A3: Correlated-disorder parameter scan"
echo ""

# Kappa scan at fixed pi=0.9
for KAPPA in 0.0 0.3 0.7 0.9; do
    for EPS in 3 6; do
        OUTFILE="${OUTDIR}/results_kappa${KAPPA}_eps${EPS}_correlated.csv"
        echo "Running: kappa=$KAPPA, eps=$EPS -> $OUTFILE"
        python3 imp_aging.py \
            --N 30 \
            --ensemble correlated \
            --epsilon $EPS \
            --n_disorder $S \
            --n_traj $M \
            --pre_sweeps $PRE \
            --beta0 $BETA0 \
            --beta $BETA \
            --tw 0 100 300 1000 3000 10000 \
            --t_max 10000 \
            --n_lags 28 \
            --kappa $KAPPA \
            --pi 0.9 \
            --seed $SEED \
            --out "$OUTFILE"
    done
done

# Pi scan at fixed kappa=0.7
for PI in 0.6 0.8 0.9 0.97; do
    for EPS in 3 6; do
        OUTFILE="${OUTDIR}/results_pi${PI}_eps${EPS}_correlated.csv"
        echo "Running: pi=$PI, eps=$EPS -> $OUTFILE"
        python3 imp_aging.py \
            --N 30 \
            --ensemble correlated \
            --epsilon $EPS \
            --n_disorder $S \
            --n_traj $M \
            --pre_sweeps $PRE \
            --beta0 $BETA0 \
            --beta $BETA \
            --tw 0 100 300 1000 3000 10000 \
            --t_max 10000 \
            --n_lags 28 \
            --kappa 0.7 \
            --pi $PI \
            --seed $SEED \
            --out "$OUTFILE"
    done
done

echo ""
echo "A3 complete!"
echo ""

# ==============================================================================
# B2: Lag-grid sensitivity (dense grid with n_lags=56)
# ==============================================================================
echo ">>> B2: Lag-grid sensitivity (n_lags=56 vs 28)"
echo ""

for EPS in 3 6; do
    for ENS in iid correlated; do
        OUTFILE="${OUTDIR}/results_dense_eps${EPS}_${ENS}.csv"
        echo "Running dense grid: eps=$EPS, ens=$ENS -> $OUTFILE"
        python3 imp_aging.py \
            --N 30 \
            --ensemble $ENS \
            --epsilon $EPS \
            --n_disorder $S \
            --n_traj $M \
            --pre_sweeps $PRE \
            --beta0 $BETA0 \
            --beta $BETA \
            --tw 0 1000 3000 10000 \
            --t_max 10000 \
            --n_lags 56 \
            --kappa 0.7 \
            --pi 0.9 \
            --seed $SEED \
            --out "$OUTFILE"
    done
done

echo ""
echo "B2 complete!"
echo ""

echo "=============================================="
echo "All robustness simulations complete!"
echo "=============================================="
echo ""
echo "Output files are in: $OUTDIR"
echo ""
echo "Next steps:"
echo "1. Run: python3 analyze_robustness_sims.py"
echo "   to generate all robustness figures"
