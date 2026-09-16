# Phase boundary: crossover incompatibility versus sequence correlation (campaigns `kappa*_T1p5`)

Five campaigns of 44 runs: `kappa` in {0, 0.25, 0.5, 0.75, 1} with `pi = 0.99`, the
11-point `eps_AB` grid {1.0, 0.975, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.3, 0.1}, four
seeds, 144 x 40 beads, `T* = 1.5` (chosen because the melt is ergodic there, so the
fluctuation estimators are meaningful), 1,000,000 production steps, box modes every
400 steps; analysis window the trailing half (2,500 tau). Code commit `b35e4be`, 8 x H100.

`delta_eps = 1 - eps_AB`; `chi = alpha delta_eps / T*` with `alpha ~ 2.0` from Stage 2.
This is the simulation's version of the (lambda, chi) phase diagram of the 2017 random
copolymer paper, with `kappa` playing the role of the sequence correlation `lambda`.

## Result 1: sequence correlation programs the transition

`S_psi(q*)` versus `delta_eps` (`S_peak_vs_delta_eps_by_kappa.png`):

| kappa | S_psi(q*) at chi = 0 | S_psi(q*) at delta_eps = 0.9 | crossover delta_eps (half-rise) | RPA spinodal (alpha = 2) |
|---|---|---|---|---|
| 0 | 1.1 | 3.2 | (no transition; 0.26 formally) | 1.50 |
| 0.25 | 2.4 | 13 | (weak; 0.27) | 0.60 |
| 0.5 | 6.2 | 121 | 0.199 | 0.216 |
| 0.75 | 15 | 502 | 0.152 | 0.104 |
| 1.0 | 25 | 1020 | 0.115 | 0.060 |

An uncorrelated random copolymer (`kappa = 0`) never demixes over the whole
incompatibility range: `S_psi(q*)` stays within a factor 3 of the ideal value 1, the
amplitudes remain Gaussian (non-Gaussianity ratio 1.7-1.9) and every mode relaxes in a
few hundred tau. `kappa = 0.25` reaches only 13. From `kappa = 0.5` upward a genuine
crossover appears, moves to weaker incompatibility as the correlation grows
(`delta_eps_c` = 0.20, 0.15, 0.115), and its plateau amplitude spans three decades
(120, 500, 1020). The mean-field spinodal reproduces the trend and agrees with the
simulation at `kappa = 0.5`, while at `kappa = 0.75-1` it lies at roughly half the
simulated crossover, i.e. mean-field theory predicts instability too early for strongly
correlated sequences. For `kappa <= 0.25` it places the spinodal at or beyond the end of
the range, consistent with no transition being seen. See `phase_boundary_kappa.png`
and `phase_boundary_kappa.csv`.

## Result 2: the second (structured) transition appears at kappa = 1

The shell anisotropy at `q*` (`max_m S_m / mean_m S_m`; 1 isotropic, 3 for a single
lamellar direction in the six-mode lowest shell) stays at 1.4-2.0 for `kappa <= 0.5`
at every incompatibility, rises to 2.2-2.4 for `kappa = 0.75` at `delta_eps >= 0.3`, and
reaches 2.7 for `kappa = 1` at `delta_eps >= 0.7`, with the peak locked at `q_min`
(`anisotropy_vs_delta_eps_by_kappa.png`). Of the two fluctuation spikes in the 2017
phase diagram (homogeneous -> random microphase, random microphase -> structured), the
first is present for `kappa >= 0.5` and the second only for the most correlated
sequences at strong incompatibility, in this box. The archived `kappa = 1`,
`eps_AB = 0.1`, `T* = 0.7` run (anisotropy 2.99 after a 1e8-step anneal) shows the
fully developed lamella.

## Result 3: peak wavevector and dynamics

For `kappa >= 0.5` the peak sits at `q* = 0.29-0.40 / sigma` and drops to `q_min` as
soon as `delta_eps >= 0.05-0.15`, i.e. at `T* = 1.5` the pattern reaches the box scale,
unlike the arrested `T* = 0.7` runs where it stays at 0.4-0.5. Beyond the crossover the
composition modes stop decaying within the 2,500-tau window (plateaus 0.9-0.97) even at
this ergodic temperature: the demixed pattern is long-lived while the melt itself is
fluid. At chi = 0 the peak relaxation time grows with `kappa` (580, 1200 tau at
`kappa = 0.5, 0.75`) because more correlated sequences carry larger, slower composition
fluctuations even without any incompatibility.

## Files

| file | contents |
|---|---|
| `phase_boundary_kappa.csv` | per kappa: crossover estimators, RPA spinodal, chi = 0 and maximum peak, anisotropy maximum |
| `phase_boundary_kappa.*` | crossover `delta_eps` versus `kappa` with the RPA spinodal |
| `S_peak_vs_delta_eps_by_kappa.*`, `anisotropy_vs_delta_eps_by_kappa.*`, `non_gaussian_ratio_vs_delta_eps_by_kappa.*` | curves per kappa |
| `kappa*_T1p5/` | per-campaign tables, `transition_summary.json`, `manifest.json`, figures |
