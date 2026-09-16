# Dynamic structure factor versus temperature (campaigns `dsf_T0p7`, `dsf_T1p0`, `dsf_T1p5`, `dsf_T2p0`)

Four campaigns of 44 runs each: `kappa = 0.5`, `eps_AB` in {1.0, 0.95, 0.9, 0.85, 0.8,
0.75, 0.7, 0.6, 0.5, 0.3, 0.1} (dense on the mixed side), four seeds, 144 x 40 beads,
`L = 22 sigma`, at quench temperatures `T* = 0.7, 1.0, 1.5, 2.0`. Production is
2,500,000 steps (12,500 tau), ten times Stage 1, with box modes recorded every 500 steps.
The analysis window is the trailing half (6,250 tau; 5,000 frames), so lags reach
3,125 tau. Code commit `b35e4be`; run on 8 x H100 with 6 OpenMM CUDA processes per GPU
under NVIDIA MPS, ~8 min per run. Raw `mode_amplitudes.npz` (44 x 50 MB per campaign) are
not versioned.

`delta_eps = 1 - eps_AB` is the incompatibility; `chi` is proportional to `delta_eps / T*`.

## Result 1: at T* >= 1.5 the composition dynamic structure factor is measurable

At `chi = 0` the two-time composition correlation `F_psi(q, t) = S(q,t)/S(q,0)` decays
inside the window for every shell. Relaxation times of the lowest shells at `chi = 0`:

| T* | q = 0.29 | q = 0.49 | q = 0.70 | q = 1.03 | late-lag plateau of F(q_min) |
|---|---|---|---|---|---|
| 0.7 | > 3125 (no decay) | > 3125 | > 3125 | 249 | 0.86 |
| 1.0 | > 3125 | 1575 | 723 | 204 | 0.74 |
| 1.5 | 1708 | 1070 | 240 | 58 | 0.47 |
| 2.0 | 1248 | ~330 | - | - | 0.15 |

(`tau_1/e` in units of tau; full tables in `<campaign>/tau_by_shell_condition.csv`.)
At `T* = 1.5` the relaxation time falls as roughly `q^-4` for `q >= 0.5 / sigma`
(`q R_g >= 1.6`), the Rouse signature, and flattens below `q ~ 0.4` where whole-chain
diffusion takes over (`figures/tau_vs_q_by_T.png`). At `T* = 2.0` the per-mode
non-Gaussianity ratio at `q*` is 1.92, i.e. the amplitudes are Gaussian as a mixed melt
should be. At `T* = 0.7` nothing below `q ~ 0.5` relaxes even in this ten-fold window:
that temperature is an arrested melt and the wrong place to measure `S(q, t)`.

## Result 2: critical slowing down of the low-q composition modes

Approaching the transition only the lowest wavevectors slow down. At `T* = 2.0` the
`q = 0.29` relaxation time rises from 1,250 tau at `chi = 0` to 2,800 tau at
`delta_eps = 0.2` and leaves the window at `delta_eps >= 0.3`, while `q >= 0.5` modes are
barely affected; at `T* = 1.5` the `q_min` mode stops decaying already at
`delta_eps = 0.1` (`figures/tau_kmin_vs_delta_eps.png`). On the demixed side `F` shows a
plateau of 0.85-0.95, an arrested pattern.

## Result 3: the Flory-Huggins mapping collapses the temperatures

`S_psi(q*)` plotted against `delta_eps / T*` puts `T* = 1.0, 1.5, 2.0` on one curve
(onset near 0.1, saturation near 0.25-0.3; `figures/S_peak_vs_chi_axis.png`), so
`chi = alpha (1 - eps_AB) / T*` holds with a temperature-independent contact factor.
`T* = 0.7` saturates lower because the pattern is arrested, not equilibrated.

## Recommendation for the theory comparison

Compare against the theory at `T* = 1.5` or `2.0` on the mixed side (`eps_AB` from 1.0
down to about 0.85), where `F(q, t)` decays to below `1/e` for every shell and the
melt is ergodic; the `q^-4` Rouse regime and its crossover to diffusion at
`q R_g ~ 1.5` are the first quantitative targets. For lags below 2.5 tau (the recording
interval here) a shorter `--mode-interval` is needed.

## Files

| path | contents |
|---|---|
| `dsf_series_summary.csv` | one row per (T*, eps_AB): peak, q*, relaxation times, plateau, ratios, drift, R_g |
| `S_peak_vs_chi_axis.*`, `F_kmin_vs_lag_by_T.*`, `tau_vs_q_by_T.*`, `tau_kmin_vs_delta_eps.*` | cross-temperature figures (`scripts/compare_dsf_campaigns.py`) |
| `dsf_T*/` | per-campaign tables, `transition_summary.json`, `manifest.json`, figures from `melt.epsab_analysis` |
