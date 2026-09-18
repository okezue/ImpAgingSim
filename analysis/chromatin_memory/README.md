# Chromatin-inspired mark memory, turnover and local feedback (campaign `chromatin_memory`)

160 runs: mark turnover `k_off` in {0.001, 0.003, 0.01, 0.03} per tau, feedback gain
`k_fb` in {0, 0.01, 0.03, 0.1, 0.3} per tau, B-B attraction `eps_BB` in {1.0, 1.5}, four
seeds. `kappa = 0.5`, initial marked (B) fraction 0.3, 144 x 40 beads at bead density 0.2
(`L = 30.65 sigma`), `eps_AA = eps_AB = 0`, `T* = 1`, 2,500,000 production steps
(12,500 tau), marks updated every 400 steps (2 tau), condensation summary every 5,000
steps, box modes every 500 steps. The basal writing rate is `k_on = k_off f_B/f_A`, so
without feedback the B fraction remains approximately 0.3. Analysis window: trailing half (6,250 tau).
Code commit `85ac983`, 8 x H100, ~10 min per run.

A and B represent unmarked and marked chromatin beads in a minimal copolymer model.
Reader-mediated attraction and writer recruitment are represented by effective B-B
attraction and the local conversion rate `k_on + k_fb H(n_B)`. Reader and writer
proteins are not explicit particles. The two attraction strengths sample the broad
aggregation crossover of the quenched-sequence campaign at `kappa = 0.5`.

[Summary figure](../../docs/figures/chromatin-memory-results.svg).

## Result 1: without feedback, B-density relaxation speeds up with turnover

At `k_fb = 0` the largest cluster holds only 2-4% of the marked beads and the dense
fraction is 1-12%. The two-time correlation of the B density at its peak wavevector,
`F_BB(q*, t)`, gives the following `1/e` relaxation times:

| k_off (1/tau) | mark memory time (tau) | B-density relaxation time, eps_BB = 1.0 (tau) | eps_BB = 1.5 (tau) |
|---|---|---|---|
| 0.001 | 477 | 93 | 290 |
| 0.003 | 205 | 51 | 94 |
| 0.01 | 67 | 41 | 48 |
| 0.03 | 23 | 16 | 21 |

The mark memory time is the measured `1/e` time of the site-mark autocorrelation.
The independent continuous-time switching prediction, `1/(k_on + k_off)`, gives
700, 233, 70 and 23 tau, respectively. The measurements approach that prediction
at faster turnover but are shorter at the slowest rate. The B-density relaxation
time also decreases with increasing turnover, and is longer at stronger attraction.
It measures collective density decorrelation, which includes mark switching and
polymer motion, rather than the tracked lifetime of an individual cluster.

## Result 2: strong feedback produces a dominant marked cluster

The largest-cluster fraction (`figures/heat_largest_cluster.png`) jumps from < 0.15 to
about 0.94-1 along a diagonal in the (`k_off`, `k_fb`) plane:

| k_off | k_fb = 0.01 | 0.03 | 0.1 | 0.3 |
|---|---|---|---|---|
| 0.001 | **0.94 / 0.99** | 1.00 | 1.00 | 1.00 |
| 0.003 | 0.06 / 0.13 | **0.95 / 0.99** | 1.00 | 1.00 |
| 0.01 | 0.02 / 0.04 | 0.04 / 0.10 | **0.97 / 1.00** | 1.00 |
| 0.03 | 0.02 / 0.03 | 0.02 / 0.03 | 0.04 / 0.11 | **0.96 / 1.00** |

(cells: `eps_BB = 1.0 / 1.5`; entries shown once round to the same value at both
attractions). On this sampled grid, `k_fb / k_off >= 10` produces a cluster containing
at least about 94% of marked beads at both attractions. Near ratio 3, the largest
cluster remains small, but B-density relaxation can take hundreds of tau. The
organization of this crossover by the rate ratio is consistent with competition
between local feedback, `k_fb H(n_B)`, and mark removal, `k_off`; its precise location
and generality require a finer scan and other model parameters.

## Result 3: feedback also increases the global marked fraction

With strong feedback (`figures/timeseries_eps_BB=1.5_fA=0.7.png`),
the B fraction rises from 0.3 to a plateau within 100-2,000 tau and the largest cluster
then contains nearly all marked beads. The plateau is governed mainly by the ratio
`k_fb / k_off`: 0.72-0.81 at ratio 10, 0.93-0.94 near ratio 30, 0.98 at ratio 100 and
0.99 at ratio 300. Attraction strength and the absolute rates still affect the
measured fractions, particularly near the crossover. Ratios of 100 or more therefore
produce near-global marking. The sampled feedback rule does not maintain a localized
marked domain at a small marked fraction. Finite writer pools, competing marks or
sequence-encoded boundaries are possible extensions for testing domain-size control.

## Result 4: density dynamics add information beyond equal-time clusters

Conditions with similarly small largest-cluster fractions can have different
`F_BB(q*, t)` relaxation times: 16-290 tau without feedback, and hundreds of tau
at some intermediate-feedback conditions. At strong feedback the `1/e` decay is
unresolved within the observation window. Jointly reporting cluster statistics,
mark autocorrelation and B-density relaxation distinguishes aggregation, mark
turnover and persistence of collective density patterns. Long-lived density
correlations alone do not demonstrate preservation of a particular spatial mark
pattern through turnover.

## Files

`per_run.csv`, `per_condition.csv` (mean and SEM over seeds of B fraction, cluster
statistics, dense fraction, `R_g`, `S_BB(q*)`, B-density relaxation time and plateau, mark memory
time, persistence fraction, runaway fraction), `summary.json`, `manifest.json`,
`figures/` (heatmaps over `k_off` x `k_fb` per `eps_BB`; B-fraction and largest-cluster
time series).
