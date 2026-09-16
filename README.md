# ImpAgingSim

Sequence correlations as a design axis for microphase structure in A/B heteropolymer melts.

The project asks a single question. If you hold composition, chemistry, chain
length and density fixed, and change only how **correlated the A/B pattern is
along the backbone**, how much does the microphase structure of the melt
change? The answer is that it changes by roughly two orders of magnitude in
the composition structure factor peak, and that a one chain sequence statistic
predicts most of it.

All simulation data is archived on Zenodo at
**[10.5281/zenodo.20499120](https://doi.org/10.5281/zenodo.20499120)**
(concept DOI, always resolves to the newest version). See
[Data availability](#data-availability).

---

## 1. The model

A melt of `M` chains, each `N` beads, in a cubic periodic box. Every bead
carries a binary label A or B. Langevin dynamics in OpenMM on GPU.

```
        periodic box, L = 22 sigma
      +-------------------------------+
      |   A-A-A-B-B-A-A      B-B-A-A  |     144 chains
      |  /           \      /       \ |     40 beads each
      | A   B-B-B-A-A       A-A-B-B   |     5,760 beads total
      |  \ /         \     /          |     rho = 0.541 sigma^-3
      |   B           B-B-B           |     T* = 0.7
      +-------------------------------+
```

### 1.1 The sequence construction

Two parameters generate the sequence. A Markov backbone sets the **range** of
correlation, and a Bernoulli mask sets its **amplitude**.

```
  step 1   Markov backbone z, persistence pi        A A A B B A A A A B B B
           eigenvalue lambda = 2*pi - 1             (long same-type runs)

  step 2   Bernoulli mask b, mean kappa             1 1 0 1 1 1 0 0 1 1 0 1
           b=1 keep backbone, b=0 redraw            . . ^ . . . ^ ^ . . ^ .

  step 3   redraw u, iid at the same f_A            . . B . . . B A . . A .

           expressed sequence                       A A B B B A A A A B A B
           psi_i = b_i z_i + (1 - b_i) u_i
```

Consequences of this construction:

- `kappa = 0` gives an independent random copolymer.
- `kappa = 1` gives the fully correlated Markov sequence.
- The lag zero variance stays independent of kappa, so the random sequence
  limit is preserved.
- The connected sequence autocorrelation scales as `kappa^2 * lambda^ell`.
- Mean composition `f_A`, chemistry, density and chain length are all held
  fixed while kappa varies. Only the covariance between labels along the chain
  changes.

Sequences are drawn independently for each chain. A Markov state never
propagates across a chain boundary.

### 1.2 The interaction potential

The nonbonded potential is split so that excluded volume is universal and
incompatibility is pair specific. This matters, because a single Lennard-Jones
term with a reduced `eps_AB` would weaken the A/B attraction and remove A/B
excluded volume at the same time, which confounds the interaction scan.

```
  U(r)
    |
    |\                  shared WCA repulsive core, identical for AA / AB / BB
    | \                 U = 4 eps_core [ (s/r)^12 - (s/r)^6 ] + eps_core
    |  \                                        for r < 2^(1/6) s
  0 +---\---------------------------------------------------- r
    |    \        ______------
    |     \______/                pair specific attractive tail
    |       |                     U = 4 eps_ab [ (s/r)^12 - (s/r)^6 ] - shift
    |   2^(1/6) s                            for r >= 2^(1/6) s
    |                                        cut and shifted at 2.5 s
```

| pair | eps | meaning |
|---|---|---|
| A-A | 1.0 | strong attraction |
| B-B | 1.0 | strong attraction |
| A-B | 0.1 | weak attraction, drives demixing |

Like attracts like, unlike barely attracts. Chain connectivity frustrates full
macrophase separation, so the melt forms microphases at a finite length scale.

The two pieces live in separate OpenMM force groups, implemented in
`melt/integrator.py`.

### 1.3 Simulation parameters

| quantity | value |
|---|---|
| chains x beads | 144 x 40 = 5,760 |
| box length | 22 sigma |
| number density | 0.541 sigma^-3 |
| quench temperature | T* = k_B T / eps = 0.7, mapped to 84.19 K |
| equilibration temperature | T* = 5.0 |
| integrator | Langevin middle, friction 1/tau, dt = 0.005 tau |
| bonds | harmonic, equilibrium length sigma |
| LJ cutoff | 2.5 sigma |
| production | 250,000 steps, snapshot every 2,000 |

Temperatures are set in kelvin using `k_B = 0.0083144626 kJ/mol/K` so that the
reduced temperature is exactly what the model intends. A regression test
asserts that T* = 0.7 maps to 84.187 K.

Initial configurations are independent random walks, then repaired by
`melt.box.relax_overlaps` before minimization. Without that repair, beads can
land about 0.007 sigma apart, and the r^-12 core puts roughly 1e25 into a
single contact, which the energy minimizer cannot always recover from. Six of
thirty runs in an early fixed density campaign diverged this way before the
repair existed.

---

## 2. What is measured

### 2.1 Static structure

The composition structure factor, taken at equal time:

```
  S_psipsi(k) = S_AA(k) + S_BB(k) - 2 S_AB(k)
```

This combination is invariant under relabelling A and B. The peak defines the
contrast `C` and the domain size `xi = 2 pi / k*`.

```
   real space, kappa = 0            real space, kappa = 1
  +-------------------+           +-------------------+
  | A B A B B A B A B |           | A A A A B B B B B |
  | B A B A A B A B A |           | A A A A B B B B B |
  | A B B A B A B A B |           | A A A A B B B B B |
  +-------------------+           +-------------------+
   interleaved, weak peak          coherent domains, strong peak

      S(k)                              S(k)
       |                                 |      /\
       |    __                           |     /  \
       |___/  \___                       |    /    \___
       +----------- k                    +--/---------- k
            k*                              k*  (lower k, larger xi)
```

Away from `f_A = 1/2` this combination mixes total density and composition
modes. The density orthogonal quantity is the Bhatia-Thornton mode
`S_cc = (1-f_A)^2 S_AA + f_A^2 S_BB - 2 f_A (1-f_A) S_AB`. Off stoichiometric
conclusions were withdrawn for this reason.

### 2.2 Dynamics

Implemented in `melt/dynamics.py`, computed from stored per bead trajectories
with periodic unwrapping:

| observable | definition |
|---|---|
| `F_s(k, t)` | self intermediate scattering function |
| `tau_alpha` | first interpolated time where `F_s = 1/e` |
| `Q(t)` | self overlap, fraction of beads within `a = 0.3 sigma` of their start |
| `chi_4(t)` | `N * Var_seeds[Q(t)]`, N being the bead count |
| `MSD` | mean square displacement |
| `alpha_2` | non Gaussian parameter |

The `chi_4` prefactor is the system size in beads. An earlier version used the
number of realizations, which underestimated the amplitude by a factor of
`N / n_runs`, about 1920 for the production geometry. A regression test pins
the normalization.

---

## 3. Results

### 3.1 Sequence correlation amplifies contrast

Raising kappa from 0 to 1 at fixed chemistry and composition:

| control | amplification of peak contrast |
|---|---|
| baseline | 18.7x |
| dense (L = 17 sigma) | 11.9x |
| soft (all attractions x 0.4) | 15.0x |
| short chain (N = 12) | 15.0x |

The domain length grows by a common factor of 2.33 across all four scans, and
the potential energy falls with kappa, consistent with more like-like contacts.

### 3.2 A one chain statistic predicts it

The finite chain intramolecular composition form factor for a Gaussian
reference chain with segment length `b`:

```
  S_0(k) / [4 f_A (1 - f_A)]  =  1 + 2 kappa^2 SUM_{ell=1}^{N-1}
                                  (1 - ell/N) lambda^ell exp(-k^2 b^2 ell / 6)
```

Plotting excess contrast against this predictor collapses the data with
logarithmic slope 1.02 and R^2 = 0.80. The melt amplifies a bare one chain
statistic. This is the same object that RPA and SCFT are built on, with
`eps_AB` playing the role of the Flory-Huggins `chi`.

### 3.3 Composition

The symmetric peak is largest at balanced composition, amplified 23.9x between
kappa = 0 and kappa = 1 at `f_A = 0.5`. Away from balance the minority species
caps the available composition variance.

### 3.4 Cross attraction

Contrast is non monotonic in `eps_AB`. The small `eps_AB` interval reads as a
broad plateau, and the earlier claim of a resolved optimum there was withdrawn. A blocked permutation test over the correlated
sequence values at `eps_AB <= 0.2` gives Q = 6.0667 and p = 0.5626, so the
interval is unresolved. Driving `eps_AB` to zero makes the melt lose cohesion
and break into globules, so it stops producing a clean microphase pattern.

### 3.5 Dynamics slow without a clean aging law

Relaxation slows with kappa and `chi_4` rises, so the correlated melt relaxes
more slowly and more heterogeneously. Relaxation times do not grow
monotonically with waiting time, and static contrast shows no consistent drift
over 1e6 to 1e8 steps, so the stored trajectories do not establish a waiting
time dependent glass transition.

Measured at each condition's own peak wavevector, the kappa = 1 to kappa = 0
ratio of `tau_alpha` is about 6.2. At a common `k = 1/sigma` it falls to about
1.9. Part of the apparent slowdown therefore reflects the larger length scale being
probed, so the structural relaxation slows more than the local bead motion does.

### 3.6 Finite size

A fixed density series at `M = 144 / 288 / 576` chains, with
`L(M) = 22 (M/144)^(1/3)` so bead density and mesh spacing stay constant,
five seeds per condition, 30 runs total:

| M | L (sigma) | kappa = 0 | kappa = 1 |
|---|---|---|---|
| 144 | 22.00 | 6.23 +/- 1.22 | 541 +/- 13 |
| 288 | 27.72 | 6.60 +/- 0.75 | 773 +/- 59 |
| 576 | 34.92 | 6.83 +/- 1.25 | 575 +/- 40 |

The roughly 90x amplification is present at every system size and `k*` is
stable. At kappa = 0 the amplitude is flat to 9%. At kappa = 1 the M = 144 and
M = 576 values differ by 0.8 sigma, so there is no evidence of systematic
drift with box size.

### 3.7 What is deliberately not claimed

The data does not establish an order-disorder transition, a universal density
exponent, or a glassy aging transition. The finite box admits a discrete set of
wavevectors, the density scans mix overlap and scattering amplitude changes,
and locating an ODT would require a temperature or `eps_AB` scan with finite
size scaling and Binder cumulants across several box sizes.

---

## 4. Repository layout

```
melt/                       the simulation engine
  integrator.py             OpenMM system: bonds, WCA core, pair specific tail
  sequences.py              Markov backbone plus Bernoulli mask generator
  box.py                    chain placement and overlap repair
  observables.py            structure factors and peak extraction
  density.py                density fields on a grid
  density_slice.py          high resolution real space slices
  direct_structure.py       direct reciprocal shell estimator
  dynamics.py               F_s, MSD, alpha_2, Q, chi_4, tau_alpha
  twotime.py                two time correlation helpers
  run.py                    single run driver
  scan.py  kappa_scan.py  temperature_scan.py  big_run.py
  fixed_density_size_scan.py    fixed density campaign driver
  fixed_density_analysis.py     deterministic aggregation for that campaign
  analyze.py  deep_analysis.py  viz.py  io.py  model.py

analysis/  analysis_aws/    derived tables and figures
aws/                        EC2 campaign scripts, see aws/README.md
tests/                      33 melt tests, 30 fixed density tests
output/                     run outputs, large directories are gitignored
archive/single_chain_mc/    superseded code, see below
```

### archive/single_chain_mc/

The original single chain Metropolis Monte Carlo study of the IMP
heteropolymer, using pair specific quenched disorder `eta_ij` and contact
overlap observables. It preceded the multi chain off lattice melt model and is
superseded by `melt/`. The code is kept so the earlier numbers stay
reproducible, and its outputs live in `output/results.csv`,
`output/results_extended.csv` and `output/robustness/`. New work belongs in
`melt/`.

The Obsidian vault in `obsidian-notes/` is local only and is gitignored.

---

## 5. Running it

### Single run

```bash
python3 -m melt.run \
  --out output/melt/demo --run_id k1.0 \
  --sequence correlated --kappa 1.0 --pi 0.99 --f_A 0.5 \
  --n_chains 144 --chain_length 40 --box_size 22.0 \
  --T_equilibrate 5.0 --T_quench 0.7 --lj_eps_AB 0.1 \
  --equilibration 30000 --n_steps 250000 --snapshot_interval 2000 \
  --compute_density --grid_size 56 --save_trajectory \
  --seed 1 --platform CUDA
```

Runs are idempotent. Re-invoking skips any run whose `structure_factor.npz`
already exists.

### Fixed density campaign

```bash
python3 -m melt.fixed_density_size_scan --dry-run              # print the 30 run plan
python3 -m melt.fixed_density_size_scan --manifest-only --platform CUDA
python3 -m melt.fixed_density_size_scan --run --platform CUDA
python3 -m melt.fixed_density_size_scan --analyze --platform CUDA
```

Each run is staged, hashed, then atomically promoted, so re-running verifies
hashes and skips completed work. A campaign is pinned to its creation commit
and requires a clean git tree. For a scheduler array, create the manifest once
and pass `--run-index "$SLURM_ARRAY_TASK_ID"`.

### Tests

```bash
python3 -m pytest tests/ -v
```

The suite covers reduced unit temperature conversion, force group separation in
the nonbonded potential, independent per chain sequence generation under the
unified `lambda = 2 pi - 1` generator, symmetric composition structure factors,
and the `chi_4 = N Var(Q)` normalization.

### AWS

Campaigns run on EC2 GPU instances. See [`aws/README.md`](aws/README.md) for
provisioning, campaign kinds, cost and teardown. Dispatch is by `SCAN_KIND`:

```bash
S3_BUCKET=your-bucket INSTANCE_TYPE=g5.2xlarge \
GIT_REF=master SCAN_KIND=fixed_density ./aws/launch.sh
```

Available kinds: `smoke`, `kappa`, `temperature`, `big`, `rerun`,
`seed_extension`, `fixed_density`, `all`.

---

## Data availability

The repository holds the code and small derived tables. The raw simulation
output is too large for git and lives on Zenodo.

**DOI [10.5281/zenodo.20499120](https://doi.org/10.5281/zenodo.20499120)**
(concept DOI, resolves to the newest version, currently 3.0.0, 3.5 GB)

| file | size | contents |
|---|---|---|
| `heteropolymer_microphase_data.tar` | 3,328 MB | full raw output of the production campaigns, mirroring the AWS results bucket. Per run `meta.json`, `snapshots.csv`, `structure_factor.npz`, and `trajectory.npz` where applicable |
| `complete_local_archive.tar` | 106 MB | the earlier single chain study, robustness sweeps, development runs, manuscript builds and working notes |
| `fixed_density_campaign.tar` | 37 MB | the 30 run fixed density finite size study, plus the superseded first execution kept for provenance |
| `Supplementary_Data_1_source_tables.zip` | 25 MB | per figure source and sensitivity tables, plus figure generation scripts |

A fresh clone gives you the code and the derived tables. Reproducing figures
from raw trajectories requires pulling the tarballs from the DOI.
