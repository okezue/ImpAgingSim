# Stage 2 eps_AB refinement at kappa = 0.5, T* = 0.7 (campaign `epsab_stage2_kappa05`)

400 runs: `eps_AB` from 0.95 to 0.65 in steps of 0.0125 (`delta_eps = 1 - eps_AB` from 0.05
to 0.35), eight seeds, two boxes at fixed bead density (144 chains, `L = 22`; 288 chains,
`L = 27.7`), 1,000,000 production steps (four times Stage 1), box modes every 400 steps.
Analysis window: trailing 500,000 steps (2,500 tau). Code commit `b35e4be`; run on two
8 x H100 nodes (~4 min per 144-chain run, ~14 min per 288-chain run). Raw output is not
versioned.

## Result 1: the mixed -> demixed crossover is at eps_AB = 0.85 +/- 0.03

`S_psi(q*)` rises from 16-19 at `delta_eps = 0.05` to a plateau of 52-55 for
`delta_eps >= 0.25`, and the half-rise point is `delta_eps = 0.148` for 144 chains and
`0.143` for 288 chains (`eps_AB = 0.852` and `0.857`). The two box sizes agree within
their standard errors at every `delta_eps` (`figures/S_peak_vs_delta_eps.png`), so the
crossover is not a finite-size artifact.

| estimator | M = 144 | M = 288 |
|---|---|---|
| half-rise of `S_psi(q*)` | 0.148 | 0.143 |
| steepest rise of `ln S_psi(q*)` (bootstrap 16-84%) | 0.169 (0.08-0.26) | 0.144 (0.144) |
| maximum of peak-intensity time variance | 0.125 (0.125-0.24) | 0.175 (0.11-0.20) |
| steepest rise of coarse-grained variance | 0.069 (0.07-0.13) | 0.094 (0.07-0.12) |
| RPA spinodal from the four most-mixed points | 0.095 (alpha = 2.09) | 0.095 (alpha = 1.96) |

The mean-field (RPA) spinodal sits at `delta_eps = 0.095` (`eps_AB = 0.905`) with a
contact factor `alpha ~ 2.0` that is now consistent between box sizes (Stage 1 with four
seeds gave `alpha` between 0.7 and 1.9 depending on wavevector). The simulation crossover
lies about 50% deeper in incompatibility than the mean-field instability, the expected
direction and size of the fluctuation correction for a finite random copolymer melt.

## Result 2: the microphase has a size-independent domain spacing

The peak wavevector stays at `q* = 0.40-0.50 / sigma` for both boxes (domain spacing
`2 pi / q* ~ 13-15 sigma`) even though the 288-chain box admits `q_min = 0.227`. The
pattern is therefore a genuine random microphase with its own length scale rather than a
box-limited macrophase. It is still coarsening slowly (mean drift +10-20% between the two
halves of the window), so plateau amplitudes are lower bounds; the crossover location on
the mixed side, where relaxation is fast, is not affected.

## Result 3: no second (structured) transition at kappa = 0.5 in this range

The shell anisotropy stays at 1.9-2.5 for both boxes across all `delta_eps`, far from the
lamellar value (3 for the six-mode lowest shell) that the archived `kappa = 1` run shows,
and the archived `kappa = 0.5`, `eps_AB = 0.1` run is isotropic (1.14). Of the two
fluctuation spikes in the 2017 phase diagram, only the homogeneous -> random microphase
one is present at `kappa = 0.5`; the random microphase -> structured transition needs
stronger sequence correlation.

## Caveat on the fluctuation estimators at T* = 0.7

The peak-intensity time variance is nearly flat (0.01-0.05) and the per-mode
non-Gaussianity ratio sits at 1.1-1.25 everywhere, because at `T* = 0.7` the composition
modes do not relax within the window (plateaus 0.90-0.96). The "spike" in density
fluctuations that marks the transition in an ergodic system is visible only in the
higher-temperature series (`analysis/dsf_series/`), where the non-Gaussianity ratio
falls from 1.9 (Gaussian, mixed) through 1.5 near `delta_eps / T* ~ 0.1`.

## Files

Same layout as `analysis/epsab_stage1/`: `per_run.csv`, `per_condition.csv`,
`spectra_by_condition.csv`, `F_peak_by_condition.csv`, `tau_by_shell_condition.csv`,
`transition_summary.json`, `manifest.json`, `figures/`.
