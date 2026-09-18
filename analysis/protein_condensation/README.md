# Sequence-programmed condensation with asymmetric energetics (campaign `protein_condensation`)

360 runs: B-B attraction `eps_BB` in {0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3}, sequence
correlation `kappa` in {0, 0.25, 0.5, 0.75, 1} (`pi = 0.99`), hydrophobic (B) fraction 0.3
and 0.5 (`f_A` = 0.7, 0.5), four seeds. 144 x 40 beads at bead density 0.2
(`L = 30.65 sigma`, implicit solvent), `eps_AA = eps_AB = 0` (only B-B pairs have
an attractive tail; all pairs share the WCA core), `T* = 1`, 30,000 equilibration steps at `T* = 3`,
1,000,000 production steps; condensation summary every 2,000 steps, box modes every
400 steps. Analysis window: trailing half. Code commit `85ac983`, 8 x H100, ~4 min per run.

B represents hydrophobic segments and A polar segments in a minimal copolymer model.
Solvent and amino-acid-specific chemistry are implicit. Two B beads are in contact
within `1.5 sigma`; clusters are connected components of this B contact network.

[Summary figure](../../docs/figures/protein-condensation-results.svg).

## Result 1: sequence correlation promotes large B-rich clusters

Fraction of B beads in the largest cluster (`figures/heat_largest_cluster.png`):

| kappa | eps_BB = 0.5 | 1.0 | 2.0 | 3.0 | (B fraction 0.3) |
|---|---|---|---|---|---|
| 0 | 0.01 | 0.02 | 0.05 | 0.06 | many small clusters, reduced chain dimensions |
| 0.25 | 0.02 | 0.06 | 0.07 | 0.07 | same |
| 0.5 | 0.06 | 0.38 | 0.32 | 0.21 | tens of medium clusters |
| 0.75 | 0.47 | 0.66 | 0.40 | 0.43 | large clusters |
| 1 | 0.79 | 0.53 | 0.48 | 0.62 | a dominant cluster |

At B fraction 0.3, uncorrelated sequences have a largest-cluster fraction of only 6%
even at `eps_BB = 3`, with about 40 clusters. Their mean radius of gyration decreases
from 3.85 to 3.47 across the attraction scan, consistent with chain compaction.
For blocky sequences (`kappa >= 0.75`), the largest cluster already holds 47-79% of B
at `eps_BB = 0.5`, while `R_g` remains about 3.73-3.75. Both sequence correlation and
attraction matter: at `kappa = 0.5`, the largest-cluster fraction changes from 0.06 to
0.38 and then to 0.21 as `eps_BB` rises from 0.5 through 1 to 3.

At B fraction 0.5, highly correlated sequences yield still larger dominant clusters
(up to 98% of B). At `kappa = 0.5`, the largest-cluster fraction reaches 0.6-0.7 for
`eps_BB = 0.75-1.75` and falls at stronger attraction. These are finite-time cluster
statistics; they identify a broad change in aggregation across the sequence scan,
without determining an equilibrium phase boundary.

## Result 2: local densification and large-cluster formation differ

The dense fraction (B beads with at least six B neighbors, `figures/heat_dense_fraction.png`)
increases across the attraction scan, from near zero to 0.7-0.99. The largest-cluster
fraction has a different, nonmonotonic dependence on attraction and changes strongly
with sequence correlation. Thus dense local B environments need not belong to one
large cluster. Weakly correlated sequences show more chain compaction and retain many
clusters; strongly correlated sequences can form large clusters while retaining
larger chain dimensions. These measurements are consistent with different balances
of intrachain compaction and interchain aggregation; distinguishing their contributions
requires resolving the chain membership of individual contacts and clusters.

## Result 3: stronger attraction can reduce the largest-cluster fraction

For several correlated-sequence conditions, stronger attraction produces a smaller
largest-cluster fraction. At `kappa = 1`, B fraction 0.3, it falls from 0.79 at
`eps_BB = 0.5` to 0.48 at 2, then rises to 0.62 at 3. At `kappa = 0.5`, the fraction
falls from 0.38 at `eps_BB = 1` to 0.21 at 3. The nonmonotonic response is
consistent with hindered coarsening at strong
attraction, but the 5,000-tau runs and aggregate cluster statistics alone do not
establish kinetic arrest, fragmentation of a pre-existing condensate, or the
equilibrium morphology.

## Files

`per_run.csv`, `per_condition.csv` (mean and SEM over seeds of cluster fraction, cluster
count, coordination numbers, dense fraction, `R_g`, `S_BB(q*)`, `F_BB(q*, t)` relaxation
time and plateau), `summary.json`, `manifest.json`, `figures/` (heatmaps over `eps_BB`
x `kappa` per composition and cluster-fraction time series).
