# Epigenetic memory: transient blobs, turnover and reader-writer feedback (campaign `chromatin_memory`)

160 runs: mark turnover `k_off` in {0.001, 0.003, 0.01, 0.03} per tau, feedback gain
`k_fb` in {0, 0.01, 0.03, 0.1, 0.3} per tau, B-B attraction `eps_BB` in {1.0, 1.5}, four
seeds. `kappa = 0.5`, initial marked (B) fraction 0.3, 144 x 40 beads at bead density 0.2
(`L = 30.65 sigma`), `eps_AA = eps_AB = 0`, `T* = 1`, 2,500,000 production steps
(12,500 tau), marks updated every 400 steps (2 tau), condensation summary every 5,000
steps, box modes every 500 steps. The basal writing rate is `k_on = k_off f_B/f_A`, so
without feedback the B fraction stays at 0.3. Analysis window: trailing half (6,250 tau).
Code commit `85ac983`, 8 x H100, ~10 min per run.

The two attractions were chosen from the protein campaign to sit at and just above the
condensation boundary for `kappa = 0.5`, i.e. the system is "not in a condition to
phase separate outright": with quenched marks a few tens of small clusters form and
the largest holds a third of the B beads at most.

## Result 1: without feedback, blobs are transient and turnover sets their lifetime

At `k_fb = 0` the largest cluster holds only 2-4% of the marked beads and the dense
fraction is 1-12%: clusters nucleate and dissolve continuously. The two-time correlation
of the B density at the cluster wavevector, `F_BB(q*, t)`, decays with

| k_off (1/tau) | mark memory time (tau) | blob lifetime, eps_BB = 1.0 | eps_BB = 1.5 |
|---|---|---|---|
| 0.001 | 477 | 93 | 290 |
| 0.003 | 205 | 51 | 94 |
| 0.01 | 67 | 41 | 48 |
| 0.03 | 23 | 16 | 21 |

The mark memory time is the `1/e` time of the site-mark autocorrelation and equals the
kinetic expectation `1/(k_on + k_off)`. Blob lifetimes track it: faster turnover,
shorter-lived blobs. At the weaker attraction the blobs die before the marks forget
(lifetime well below the memory time: they dissolve by polymer diffusion); at the
stronger attraction the two times converge and mark turnover becomes the limiting step.

## Result 2: feedback stabilizes a blob once k_fb / k_off exceeds about 3-10

The largest-cluster fraction (`figures/heat_largest_cluster.png`) jumps from < 0.15 to
> 0.95 along a diagonal in the (`k_off`, `k_fb`) plane:

| k_off | k_fb = 0.01 | 0.03 | 0.1 | 0.3 |
|---|---|---|---|---|
| 0.001 | **0.94 / 0.99** | 1.00 | 1.00 | 1.00 |
| 0.003 | 0.06 / 0.13 | **0.95 / 0.99** | 1.00 | 1.00 |
| 0.01 | 0.02 / 0.04 | 0.04 / 0.10 | **0.97 / 1.00** | 1.00 |
| 0.03 | 0.02 / 0.03 | 0.02 / 0.03 | 0.04 / 0.11 | **0.96 / 1.00** |

(cells: `eps_BB = 1.0 / 1.5`). Stabilization needs `k_fb / k_off >= 10` at both
attractions; the ratio 3 is marginal (largest cluster 0.06-0.13, blob lifetime raised to
200-900 tau, above the mark memory time), which is the transition region. Below the
threshold the outcome is the transient regime of Result 1 regardless of attraction. The
condition is a ratio because a blob survives when writing inside it (`~ k_fb H(n_B)`,
with `H ~ 0.5-0.6` for the 6-8 marked neighbors in a dense region) outruns the
turnover `k_off` of its own marks.

## Result 3: a stabilized blob does not stay a blob

Once stabilized the domain accretes marks (`figures/timeseries_eps_BB=1.5_fA=0.7.png`):
the B fraction rises from 0.3 to a plateau within 100-2,000 tau and the largest cluster
then contains essentially every marked bead. The plateau depends only on the ratio
`k_fb / k_off`: 0.72-0.81 at ratio 10, 0.93-0.94 at ratio 30, 0.98 at ratio 100, 0.99 at
ratio 300, identical across `k_off` and nearly identical across `eps_BB`. For
`k_fb / k_off >= 100` the whole system converts to B (runaway). In this minimal model
positive feedback therefore has no intrinsic size control: the same condition that
rescues a transient blob makes it grow until turnover at the domain surface balances
writing inside, and that balance is set by the rate ratio, not by the polymer. A
finite, stable heterochromatin domain needs an ingredient the model does not yet
contain, such as a finite writer pool, a competing mark that stabilizes A, or
sequence-encoded boundaries. This is the main open point to take to the theory.

## Result 4: the two-time correlation is the right diagnostic

The equal-time cluster statistics of the transient regime and of the marginal regime
are similar (largest cluster 0.02-0.13), but the blob lifetime from `F_BB(q*, t)`
separates them cleanly: 16-290 tau (below or at the mark memory time) when blobs are
transient, 200-900 tau (above it) when feedback is starting to hold them, and beyond
the window when they are stabilized. As anticipated in the chromatin note, the
two-time density correlation, not the structure factor, tells a fluctuation from a
domain.

## Files

`per_run.csv`, `per_condition.csv` (mean and SEM over seeds of B fraction, cluster
statistics, dense fraction, `R_g`, `S_BB(q*)`, blob lifetime and plateau, mark memory
time, persistence fraction, runaway fraction), `summary.json`, `manifest.json`,
`figures/` (heatmaps over `k_off` x `k_fb` per `eps_BB`; B-fraction and largest-cluster
time series).
