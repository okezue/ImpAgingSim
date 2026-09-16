# Stage 1 eps_AB sweep at kappa = 0.5 (campaign `epsab_kappa05`)

148 runs: `eps_AB` from 1.0 to 0.1 in steps of 0.025, four seeds each, 144 x 40 beads,
`L = 22 sigma`, `kappa = 0.5`, `pi = 0.99`, `f_A = 1/2`, exact global composition,
`T* = 0.7`, 30,000 equilibration steps at `T* = 5`, 250,000 production steps.
Analysis window: trailing half of production (125,000 steps, 625 frames of box modes).
Code commit `66b24ff`; design hash in `manifest.json`. Raw output (`mode_amplitudes.npz`,
`snapshots.csv`, `structure_factor.npz` per run, 148 x ~15 MB) is not versioned.

The incompatibility axis is `delta_eps = eps_AA - eps_AB = 1 - eps_AB`.

## Result 1: the chi = 0 reference is mixed

At `eps_AB = 1` the composition spectrum is flat, `S_psi(q) = 3-5` for `q <= 1/sigma`,
against the ideal single-chain form factor `S_0(q_min) = 6.9`, `S_0(0.7) = 3.6` for
`b = 1.24 sigma` (from the measured `R_g = 3.21`). Entropy of mixing wins, as expected.

## Result 2: demixing sets in between delta_eps ~ 0.1 and 0.3

`S_psi(q*)` rises continuously from 6 (delta_eps = 0) through 10 (0.025-0.075),
16 (0.1), 20 (0.15-0.175), 22-29 (0.2-0.25), 30-38 (0.3-0.4) and then saturates at
40-48 for `delta_eps >= 0.4`, with a peak at `q* ~ 0.5-0.6 / sigma` (domain spacing
`2 pi / q* ~ 11 sigma`, half the box). See `figures/S_peak_vs_delta_eps.png` and
`figures/spectra_M144.png`.

| estimator | delta_eps | seed bootstrap (16-84%) | eps_AB |
|---|---|---|---|
| max of peak-intensity time variance | 0.20 | 0.075 - 0.20 | 0.80 |
| steepest rise of `ln S_psi(q*)` | 0.0125 | 0.0125 - 0.0875 | ~0.99 |
| RPA spinodal, alpha fitted on the four most-mixed points at `q = 0.45` | 0.107 | (alpha = 1.89) | 0.89 |
| RPA spinodal, six points at `q = 0.40 / 0.49 / 0.64` | 0.30 / 0.15 / 0.11 | (alpha = 0.67 / 1.36 / 1.84) | 0.70 / 0.85 / 0.89 |

Reading: fluctuations are amplified from the very first step away from chi = 0 (the
mean-field regime), the transition to a demixed pattern sits in
`0.1 <= delta_eps <= 0.3` (`0.7 <= eps_AB <= 0.9`), and the RPA spinodal falls at the
low end of that bracket, `delta_eps ~ 0.1-0.15`. Mean-field theory locating the
instability slightly before the real system demixes is the expected direction of the
fluctuation correction for random copolymers. The bridge `chi = alpha (1 - eps_AB)/T*`
has `alpha ~ 1-2`; the spread across wavevectors is dominated by seed noise on the
mixed side (four seeds, six modes in the lowest shell), which Stage 2 must fix.

## Result 3: the 250k-step protocol does not equilibrate the demixed side

The archived kappa = 0.5, `eps_AB = 0.1` runs on Zenodo that were annealed for 1e8 steps
reach `S_psi(q_min) ~ 160` with the peak at `q_min`; the same condition here after
250k steps reaches 48 with the peak at `q = 0.64`. 72% of runs show a positive drift
of `S_psi(q*)` between the two halves of the window (mean +6.5%). The plateau at
40-48 is therefore a coarsening transient limited by run length, not an equilibrium
amplitude, and the onset region is where a longer production run matters most
(critical slowing down).

## Result 4: composition dynamics are slower than the window at T* = 0.7

`F_psi(q*, t)` decays only to 0.95 by `t = 300 tau` at every `eps_AB`, including
chi = 0. Only shells with `q >= 1/sigma` reach `1/e` inside the window (at
`tau ~ 180-300 tau`). The temporal non-Gaussianity ratio (~1.1) and the shell
anisotropy (~2.4, the value expected for six unrelaxed modes) are therefore
unsampled null values here, not physics. Measuring the dynamic structure factor's
decay requires either production runs 10-100x longer with sparser mode recording, or
a higher temperature (`T* = 1-2`) where the melt is an ergodic liquid; in both cases
the chi axis stays `(1 - eps_AB)/T*`. See `figures/F_peak_vs_lag_M144.png`.

## Files

| file | contents |
|---|---|
| `per_run.csv` | one row per run: peak, spectra summaries, fluctuation estimators, relaxation fits, energy, `R_g` |
| `per_condition.csv` | mean and SEM over seeds per `eps_AB` |
| `spectra_by_condition.csv` | `S_psi(q)` per shell and `eps_AB`, mean and SEM over seeds |
| `F_peak_by_condition.csv` | `F_psi(q*, t)` per `eps_AB`, mean and SEM over seeds |
| `tau_by_shell_condition.csv` | `1/e` relaxation time per shell and `eps_AB` |
| `transition_summary.json` | transition estimators with bootstrap, RPA fit |
| `manifest.json` | campaign design, settings, commit and run list |
