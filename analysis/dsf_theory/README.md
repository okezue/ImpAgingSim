# Theory-comparison dataset: composition dynamic structure factor on the mixed side

Campaigns `dsf_theory_T1p5` and `dsf_theory_T2p0` (56 runs each): `kappa = 0.5`,
`eps_AB` in {1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85}, eight seeds, 144 x 40 beads,
`L = 22 sigma`, 5,000,000 production steps (25,000 tau) with box modes every 200 steps
(1 tau). Analysis window: trailing half (12,500 tau; 12,500 frames), lags to 6,250 tau.
`dsf_chi0_M288_T1p5` (4 runs): the chi = 0 point in the 288-chain box (`L = 27.7`,
`q_min = 0.227`) with the same protocol. Code commit `b35e4be`, 8 x H100.

This is the dataset intended for comparison with the dynamic theory: an ergodic melt,
the mixed side up to the crossover, every shell relaxing inside the window, and enough
time resolution to see the initial decay for `q <= 1.5 / sigma`.

## Result 1: Rouse scaling of the composition relaxation at chi = 0 (T* = 1.5)

`tau_1/e(q)` of `F_psi(q, tau) = S(q,tau)/S(q,0)` at `eps_AB = 1`:

| q (1/sigma) | 0.29 | 0.40 | 0.49 | 0.57 | 0.64 | 0.70 | 0.81 | 0.90 | 0.99 | 1.14 |
|---|---|---|---|---|---|---|---|---|---|---|
| tau (tau) | 3130 | 1211 | 785 | 508 | 398 | 274 | 165 | 102 | 68 | 38 |

A power-law fit for `q >= 0.55` gives `tau ~ q^-3.92`, against the Rouse prediction
`q^-4` for `q R_g > 1` (`R_g = 3.2`, so `q R_g >= 1.8` in the fit range), and the decay at
chi = 0 is a single exponential (stretched-exponential exponent `beta = 0.99` at `q*`).
Below `q ~ 0.4` the curve bends toward the whole-chain diffusion regime, which the
288-chain runs probe down to `q = 0.227`.

## Result 2: critical slowing down is confined to the lowest wavevectors

| eps_AB | delta_eps | S_psi(q_min) | tau(q = 0.29) | tau(q = 0.40) | tau(q = 0.57) | tau(q >= 0.8) | beta at q_min |
|---|---|---|---|---|---|---|---|
| 1.0 | 0 | 7.2 | 3130 | 1211 | 508 | unchanged | 0.99 |
| 0.95 | 0.05 | 13.8 | 3756 | 2142 | 654 | unchanged | 1.22 |
| 0.90 | 0.10 | 31.8 | 5938 | 2424 | 827 | unchanged | 1.08 |
| 0.875 | 0.125 | 41.3 | 5193 | 3758 | 923 | unchanged | 0.95 |
| 0.85 | 0.15 | 59.8 | > 6250 | 3347 | 759 | unchanged | 0.61 |

Between chi = 0 and the crossover the amplitude at `q_min` grows eightfold and its
relaxation time doubles and then leaves the window, while every shell with
`q >= 0.8 / sigma` keeps the chi = 0 relaxation time to within 10%. The slowing is
therefore purely collective and long-wavelength, as a dynamic RPA of the composition
field predicts (`Gamma(q) ~ q^2 / S(q)` style thermodynamic slowing). At the crossover
the decay becomes stretched (`beta = 0.61`) and a plateau of 0.7 appears: the pattern
is beginning to arrest even in this fluid melt.

## Result 3: temperature dependence and the long-wavelength limit

At `T* = 2.0` the same picture holds with `tau ~ q^-3.82` for `q >= 0.55` at chi = 0
and `S_psi(q_min)` rising only from 7.8 to 22.6 by `eps_AB = 0.85` (against 59.8 at
`T* = 1.5`), consistent with `chi = alpha delta_eps / T*`. The chi = 0 relaxation times at
`T* = 2.0` are 2.1-2.3 times shorter than at `T* = 1.5` at every wavevector (1458 vs
3130 tau at `q_min`, 581 vs 1211 at `q = 0.40`, 344 vs 785 at `q = 0.49`), whereas a
kinetic coefficient proportional to `T*` alone would give 1.33: the bead mobility in this
dense melt has an activated component that the theory's friction constant must absorb.

The 288-chain chi = 0 runs extend the wavevector range down to `q = 0.227`, where
`tau = 3175 tau`, essentially equal to the 144-chain value at `q = 0.286` (3130); between
`q = 0.23` and `0.32` the relaxation time levels off near 3,000 tau instead of continuing
the `q^-2` growth expected from whole-chain diffusion (a free-draining estimate gives
500 tau at `q = 0.227`). The 144- and 288-chain curves agree within errors at common
wavevectors, so this is not a box artifact. Why the longest-wavelength composition mode
relaxes no slower than the `q = 0.29` mode is an open question for the theory comparison.

## Result 4: the initial decay rate

The short-time rate `Gamma = -d ln F / dt` at `q_min` falls from `1.6e-3 / tau` at
chi = 0 to `2.6e-4 / tau` at `eps_AB = 0.85`, a factor 6 for an eightfold rise in
`S(q_min)`: `Gamma S(q)` is approximately constant, the fingerprint of a
thermodynamically driven slowdown with an incompatibility-independent kinetic
coefficient. This ratio is the natural first quantity to compare with the theory.

## Files

Per campaign: `per_run.csv`, `per_condition.csv`, `spectra_by_condition.csv`,
`F_peak_by_condition.csv` (`F(q*, tau)` on the 1-tau lag grid, mean and SEM over
seeds), `tau_by_shell_condition.csv`, `transition_summary.json`, `manifest.json`,
`figures/`. The full mode-resolved `F(q, tau)` for every shell and seed is in the
per-run `dynamic_structure.npz` (in `analysis/runs/` of the raw archive) and can be
recomputed from `mode_amplitudes.npz` with `python -m melt.dynamic_structure RUN_DIR`.
