# Sequence-programmed condensation with asymmetric energetics (campaign `protein_condensation`)

360 runs: B-B attraction `eps_BB` in {0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3}, sequence
correlation `kappa` in {0, 0.25, 0.5, 0.75, 1} (`pi = 0.99`), hydrophobic (B) fraction 0.3
and 0.5 (`f_A` = 0.7, 0.5), four seeds. 144 x 40 beads at bead density 0.2
(`L = 30.65 sigma`, free volume as implicit solvent), `eps_AA = eps_AB = 0` (B loves B,
A-B neutral, shared WCA core), `T* = 1`, 30,000 equilibration steps at `T* = 3`,
1,000,000 production steps; condensation summary every 2,000 steps, box modes every
400 steps. Analysis window: trailing half. Code commit `85ac983`, 8 x H100, ~4 min per run.

B is the hydrophobic (or reader-bound) species; A is polar (or unmarked). Two beads are in
contact within `1.5 sigma`; clusters are connected components of the B contact network.

## Result 1: the pattern, not the attraction, decides whether a condensate forms

Fraction of B beads in the largest cluster (`figures/heat_largest_cluster.png`):

| kappa | eps_BB = 0.5 | 1.0 | 2.0 | 3.0 | (B fraction 0.3) |
|---|---|---|---|---|---|
| 0 | 0.01 | 0.02 | 0.05 | 0.06 | many small clusters, chains collapse individually |
| 0.25 | 0.02 | 0.06 | 0.07 | 0.07 | same |
| 0.5 | 0.06 | 0.38 | 0.32 | 0.21 | tens of medium clusters |
| 0.75 | 0.47 | 0.66 | 0.40 | 0.43 | a few large condensates |
| 1 | 0.79 | 0.53 | 0.48 | 0.62 | one dominant condensate |

The boundary is a horizontal line at `kappa ~ 0.5`. Uncorrelated hydrophobic patterns
never phase separate: even at `eps_BB = 3` the largest cluster holds 6% of B while the
number of clusters stays at ~40 and the chain radius of gyration shrinks from 3.85 to
3.47, i.e. each chain collapses onto its own scattered hydrophobic beads (intramolecular
micelles). Blocky patterns (`kappa >= 0.75`) condense already at `eps_BB = 0.5`, where
the dense fraction is still low: the blocks find each other before the attraction is
strong enough to compact them, and the chains keep their open size (`R_g` 3.65-3.75).
At B fraction 0.5 the same picture holds with a larger single condensate for
`kappa >= 0.75` (up to 98% of B) and a crossover value of `kappa = 0.5` that yields one
condensate at intermediate attraction (0.6-0.7 for `eps_BB` 0.75-1.75) but many
clusters at strong attraction.

## Result 2: two condensation routes

The dense fraction (B beads with at least six B neighbors, `figures/heat_dense_fraction.png`)
grows monotonically with `eps_BB` at every `kappa`, from ~0 to 0.7-0.99, while the
largest-cluster fraction depends on `kappa` alone. Strong attraction therefore always
produces dense B environments, but only sequence correlation decides whether those
environments are shared between chains (a condensate) or private to each chain (a
collapsed globule). For `kappa <= 0.25` the collapse route dominates: `R_g` falls by
10-20% while the cluster count stays high. For `kappa >= 0.75` the condensation route
dominates: `R_g` is nearly unchanged and the cluster count drops to 2-4.

## Result 3: strong attraction fragments the condensate

For `kappa = 0.5-1` the largest-cluster fraction is largest at moderate attraction and
falls again for `eps_BB >= 2` (e.g. 0.79 -> 0.48 for `kappa = 1`, B fraction 0.3), while
the number of clusters rises from 3 to 4-10. Strong attraction arrests the coarsening
of the early clusters before they merge, the kinetic trapping familiar from the melt
coarsening study: the equilibrium state is presumably a single condensate, but the
dense clusters stop exchanging material within the 5,000-tau window.

## Files

`per_run.csv`, `per_condition.csv` (mean and SEM over seeds of cluster fraction, cluster
count, coordination numbers, dense fraction, `R_g`, `S_BB(q*)`, `F_BB(q*, t)` relaxation
time and plateau), `summary.json`, `manifest.json`, `figures/` (heatmaps over `eps_BB`
x `kappa` per composition and cluster-fraction time series).
