# IMP Heteropolymer Aging Simulation (Metropolis MC)

Complete pipeline for simulating aging dynamics in the IMP heteropolymer with
tunable disorder correlations.

## Pipeline

| Script | Purpose |
|--------|---------|
| `imp_aging.py` | Metropolis MC simulation, temperature quench, two-time observables |
| `analyze_results.py` | Disorder-averaged timescales, aging exponent mu, collapse scores |
| `extended_analysis.py` | KWW fitting, MSD/alpha2 analysis, energy/Rg tracking |
| `quantitative_analysis.py` | Power-law fits, ensemble separation metrics |
| `robustness_analysis.py` | Bootstrap CIs, multi-threshold tau, smoothed t* |
| `create_comparison_figures.py` | Overlay figures for iid vs correlated |
| `test_imp_aging.py` | 47-test suite (energy, disorder, observables, KWW) |

## Quickstart

```bash
python3 imp_aging.py \
  --out output/results.csv \
  --ensemble iid correlated \
  --epsilon 0 3 6 \
  --n_disorder 20 --n_traj 8 \
  --seed 123

python3 analyze_results.py --csv output/results.csv --outdir analysis --mono
python3 extended_analysis.py --csv output/results.csv --outdir analysis/extended
python3 -m pytest test_imp_aging.py -v
```

## Observables

| Observable | Description |
|-----------|-------------|
| Q(t_w,t) | Contact overlap (fraction of non-bonded pairs in contact at both t_w and t_w+t) |
| chi4 | N_pairs * (var Q across trajectories) - dynamical heterogeneity |
| D4 | Mean squared change in pairwise distances |
| MSD | Per-bead mean-square displacement |
| alpha2 | Non-Gaussianity parameter (3/5)*<r^4>/<r^2>^2 - 1 |
| Energy | Total potential energy at each snapshot |
| Rg | Radius of gyration |
| n_contacts | Number of non-bonded pairs in contact |

## Disorder Ensembles

- **iid**: eta_ij ~ N(0,1) independent
- **correlated**: eta_ij = kappa * sigma_i * sigma_j + sqrt(1-kappa^2) * xi_ij,
  where sigma is a Markov chain on {+1,-1} with persistence pi
