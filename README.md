# IMP Heteropolymer Aging Simulation (Metropolis MC)

This folder contains a complete, reproducible pipeline:

1. `imp_aging.py`: runs Metropolis Monte Carlo for the IMP-style heteropolymer,
   performs a temperature quench, and records two-time (aging) observables.
2. `analyze_results.py`: post-processes the CSV output to compute disorder-
   averaged timescales/exponents and writes a publication-style figure pack
   (PNG+PDF plus a multi-page PDF).

## Quickstart

Run the simulation (this can take hours depending on statistics):

```bash
python3 imp_aging.py \
  --out results.csv \
  --ensemble iid correlated \
  --epsilon 0 3 6 \
  --n_disorder 20 \
  --n_traj 8 \
  --pre_sweeps 2000 \
  --tw 0 100 300 1000 3000 10000 \
  --t_max 10000 \
  --n_lags 28 \
  --seed 123
```

Analyze and create figures:

```bash
python3 analyze_results.py \
  --csv results.csv \
  --outdir analysis \
  --mono \
  --figure_pack \
  --zip_pack
```

Outputs:
- `analysis/summary_per_tw.csv`
- `analysis/summary_condition.csv`
- `analysis/figure_pack/` (PNG+PDF per figure, `figure_pack.pdf`, `FIGURES.csv`)
- `analysis/figure_pack.zip` (if `--zip_pack`)

## Notes

- The simulation outputs one row per `(ensemble, epsilon, disorder_idx, tw, lag)`.
- `chi4` is computed as `N_pairs * ( <Q^2> - <Q>^2 )` across trajectories.
- The analysis computes `tau_alpha(tw)` as the first crossing of
  `Q_norm(tw, t) = Q(tw, t) / Q(tw, 0)` through `exp(-1)`.
