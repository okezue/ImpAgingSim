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

![A forty-bead A/B chain beside a schematic periodic multichain melt; a highlighted chain segment repeats after translation by one box length.](docs/figures/polymer-model.svg)

**Figure 1. Bead-spring geometry.** Panel A shows a chain with harmonic bonds;
panel B shows a schematic melt. Teal and orange denote A and B. The highlighted
segment repeats across a periodic boundary separated by $L=22\sigma$.
Production systems contain 144 chains of 40 beads.

### 1.1 The sequence construction

Two parameters generate the sequence. A Markov backbone sets the **range** of
correlation, and a Bernoulli mask sets its **amplitude**.

![Aligned backbone, keep/redraw mask, fresh draws, and final bead labels, beside discrete covariance curves separating correlation amplitude from range.](docs/figures/sequence-construction.svg)

**Figure 2. Sequence construction.** In panel A, $b_i=1$ keeps the Markov
backbone label and $b_i=0$ redraws independently at the same mean A fraction:
$\psi_i=b_i z_i+(1-b_i)u_i$, with $\Pr(b_i=1)=\kappa$. A redraw can retain
the original species. Panel B separates amplitude from range. For the
unconditioned generator, the normalized connected covariance is $\Gamma(0)=1$ and
$\Gamma(\ell)=\kappa^2(2\pi-1)^\ell$ for $\ell\geq1$. The plotted values
illustrate amplitude and range separately; they are analytic, not simulation data.

Consequences of this construction:

- `kappa = 0` gives an independent random copolymer.
- `kappa = 1` gives the fully correlated Markov sequence.
- The lag zero variance stays independent of kappa, so the random sequence
  limit is preserved.
- At positive lag, the normalized connected sequence autocorrelation is
  `kappa^2 * lambda^ell`, with `lambda = 2*pi - 1`.
- Mean composition `f_A`, chemistry, density and chain length are all held
  fixed while kappa varies. Only the covariance between labels along the chain
  changes.

The ordinary generator draws chains independently and fixes composition in
expectation. A Markov state never propagates across a chain boundary. The
optional exact-total generator, used by the fixed-density campaign, conditions
these draws on the global A count.

### 1.2 The interaction potential

The nonbonded potential is split so that excluded volume is universal and
incompatibility is pair specific. This matters, because a single Lennard-Jones
term with a reduced `eps_AB` would weaken the A/B attraction and remove A/B
excluded volume at the same time, which confounds the interaction scan.

![Nonbonded A-A, B-B, and A-B pairs beside the exact implemented potential: a shared repulsive core and pair-specific attractive branches with a jump at their onset.](docs/figures/interaction-potential.svg)

**Figure 3. Implemented nonbonded energy.** Panel A identifies pair types and
attraction strengths; panel B shows their energies. The WCA core is common to
all pairs. Attraction begins at $r_m=2^{1/6}\sigma$ and is shifted to zero at
$r_c=2.5\sigma$. Open and filled endpoints show the energy jump at $r_m$
caused by the current step-gated tail. Adjacent bonded beads are excluded
from both nonbonded terms.

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

$$
S_{\psi\psi}(k)=S_{AA}(k)+S_{BB}(k)-2S_{AB}(k).
$$

This combination is invariant under relabelling A and B. Its peak defines the
contrast `C`. A resolved peak wavevector gives a characteristic spacing
$\xi=2\pi/k^*$; a maximum at the lowest sampled wavevector is limited by the box.

![Measured composition spectra for random and correlated chains, with five-seed uncertainty bands, and their peak amplitudes across three fixed-density system sizes.](docs/figures/measured-structure.svg)

**Figure 4. Measured static structure.** Panel A shows 144 chains; panel B
compares three system sizes at fixed density. All use 40 beads per chain,
$f_A=0.5$, $\pi=0.99$, and $T^*=0.7$. The direct reciprocal-shell
estimator uses the total-bead-normalized channel $S_{\psi\psi}^{(N)}/2$.
Bands and bars show mean ± SEM across five independent seeds, each averaging
its final five saved configurations. At $\kappa=1$, the two smaller boxes peak
at their lowest sampled wavevector, $k_{\min}=2\pi/L$.
[Figure sources and regeneration](docs/figures/README.md).

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

$$
\frac{S_0(k)}{4f_A(1-f_A)}
=1+2\kappa^2\sum_{\ell=1}^{N-1}
\left(1-\frac{\ell}{N}\right)\lambda^\ell
\exp\!\left(-\frac{k^2b^2\ell}{6}\right).
$$

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
five seeds per condition, 30 runs total. Entries are selected peak
$S_{\psi\psi}^{(N)}/2$, mean ± SEM, as plotted in Figure 4B:

| M | L (sigma) | kappa = 0 | kappa = 1 |
|---|---|---|---|
| 144 | 22.00 | 6.23 +/- 1.22 | 541 +/- 13 |
| 288 | 27.72 | 6.60 +/- 0.75 | 773 +/- 59 |
| 576 | 34.92 | 6.83 +/- 1.25 | 575 +/- 40 |

Peak amplification ranges from 84× to 117× across these sizes. At kappa = 0,
the endpoint amplitudes differ by about 10%. At kappa = 1, the M = 144 and
M = 576 amplitudes differ by 0.83 combined SEM, with a higher intermediate
value. The correlated peak occurs at $k_{\min}=2\pi/L$ in the two smaller
boxes and at $\sqrt{2}k_{\min}$ in the largest. This series therefore does
not establish a size-independent domain length.

### 3.7 What is deliberately not claimed

The data does not establish an order-disorder transition, a universal density
exponent, or a glassy aging transition. The finite box admits a discrete set of
wavevectors, the density scans mix overlap and scattering amplitude changes,
and locating an ODT would require a temperature or `eps_AB` scan with finite
size scaling and Binder cumulants across several box sizes.

---

## 4. The incompatibility sweep and the dynamic structure factor

Two extensions take the model to the next question: at what A/B
incompatibility does the melt cross from mixed to demixed at fixed sequence
correlation, and how does the composition pattern relax in time.

### 4.1 Dictionary to the field theory

| simulation | field theory |
|---|---|
| `kappa` | the sequence correlation parameter `lambda` of random copolymer theory |
| `eps_AB` at fixed `eps_AA = eps_BB = 1` | the Flory-Huggins `chi`; the incompatibility axis is `delta_eps = 1 - eps_AB` |
| `eps_AB = 1` | `chi = 0`: every pair interacts identically, entropy of mixing wins |
| `eps_AB = 0.1` | the strongly demixed production value of sections 1-3 |

Because the WCA core is shared by all pairs, lowering `eps_AB` weakens the A-B
attraction without touching excluded volume, so the sweep is a clean
incompatibility knob. The bridge `chi = alpha (1 - eps_AB) / T*` has one
unknown, the effective contact number `alpha`, which `melt.rpa` fits from the
mixed side of the sweep where mean-field theory applies.

### 4.2 The sweep (`melt.epsab_scan`)

Everything is held at the production protocol (144 x 40 beads, `L = 22`,
`kappa = 0.5`, `pi = 0.99`, `f_A = 1/2`, exact global composition,
`T* = 0.7`, 30,000 equilibration steps at `T* = 5`, 250,000 production steps)
and only `eps_AB` moves: 37 values from 1.0 down to 0.1 in steps of 0.025,
four seeds each, 148 runs. The grid is deliberately dense because a
finite-system transition is sudden and its location is not known in advance;
a second stage refines the bracketed interval with a finer list, more seeds
and a second box size under a new campaign id.

Every run records the exact Fourier amplitudes of the A and B bead densities
at all periodic-box modes with `|q| <= 1.5 / sigma` (618 modes in 24 shells
for `L = 22`) every 200 steps, so 1,250 frames per run. The evaluation uses
the per-axis factorization of `exp(i q . r)` for integer box modes and is
exact; recording leaves the seeded trajectory bit-identical.

### 4.3 What is measured

Static, from the trailing half of production (the stationary window):

| observable | definition | reads |
|---|---|---|
| `S_psi(q*)` | time-averaged peak of `S_AA + S_BB - 2 S_AB` | amplitude of composition fluctuations |
| coarse-grained variance | `(1/N) sum_q S_psi(q) exp(-q^2 l^2)` | real-space density fluctuation at scale `l` |
| peak-intensity variance | `Var_t[S_psi(q*, t)] / <S_psi(q*)>^2` | fluctuations of the fluctuations; spikes at a transition |
| non-Gaussianity ratio | `<|rho_psi(q*)|^4> / <|rho_psi(q*)|^2>^2` | 2 for Gaussian (mixed) amplitudes, 1 for a frozen pattern |
| seed variance of `S_psi(q*)` | relative variance across seeds | disorder-to-disorder susceptibility |

Dynamic, the two-time object Spakowitz asked for:

$$
S_{\psi\psi}(q,\tau)=\frac{1}{N}\left\langle \rho_\psi(q,t_0)\,\rho_\psi^*(q,t_0+\tau)\right\rangle_{t_0,\ |q|\in\text{shell}},
\qquad F(q,\tau)=\frac{S(q,\tau)}{S(q,0)} .
$$

It is computed for every shell by the Wiener-Khinchin identity over all time
origins in the window, together with the total-density channel
`rho_A + rho_B`. Time means are not subtracted, so an arrested pattern
appears as a plateau in `F(q, tau -> inf)` rather than being removed. Per
shell the analysis reports the `1/e` relaxation time, the short-time decay
rate from `-ln F`, a stretched-exponential fit and the late-lag plateau.

### 4.4 Locating the transition

`melt.epsab_analysis` aggregates a campaign into per-run and per-condition
tables (mean and SEM over seeds), spectra and `F(q*, tau)` per condition,
and `transition_summary.json` with five estimators of the transition along
`delta_eps`, each with a seed bootstrap: the maximum of the peak-intensity
variance, the maximum of the seed variance, the steepest rise of
`ln S_psi(q*)`, the steepest rise of the coarse-grained variance, and the
crossing of the non-Gaussianity ratio through 1.5. The same file holds the
RPA comparison: `alpha`, the fitted `chi(eps_AB)`, the predicted spinodal
`eps_AB` and the mean-field `S(q*)` curve that the figures overlay on the
simulation.

The RPA used is the incompressible one-component form in the simulation's
normalization, `S^{-1} = S_0^{-1} - chi/2`, with `S_0(q)` the finite-chain
composition form factor of section 3.2 and the segment length taken from the
measured `R_g` at `eps_AB = 1`. For a symmetric blend it reproduces
`chi_s N = 2`.

---

## 5. Repository layout

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
  modes.py                  exact box-mode amplitudes rho_A(q,t), rho_B(q,t) recorder
  dynamic_structure.py      S(q,tau), F(q,tau) and static fluctuation observables
  rpa.py                    finite-chain S_0(q), RPA S(q), spinodal, eps_AB -> chi bridge
  run.py                    single run driver
  scan.py  kappa_scan.py  temperature_scan.py  big_run.py
  campaign.py               shared restartable-campaign machinery (manifest, hashes, locks)
  epsab_scan.py             eps_AB sweep at fixed kappa, shardable and parallel
  epsab_analysis.py         sweep aggregation, transition estimators, RPA comparison, figures
  fixed_density_size_scan.py    fixed density campaign driver
  fixed_density_analysis.py     deterministic aggregation for that campaign
  analyze.py  deep_analysis.py  viz.py  io.py  model.py

docs/figures/              README SVGs, PNG exports, and figure provenance
scripts/render_*_figure*   reproducible README figure generators
scripts/fetch_zenodo.py    download and verify the Zenodo archive
scripts/modes_from_trajectory.py   archived trajectory.npz -> mode_amplitudes.npz
analysis/  analysis_aws/    derived tables and figures
aws/                        EC2 campaign scripts, see aws/README.md
sherlock/                   SLURM kit for Stanford Sherlock, see sherlock/README.md
tests/                      113 tests across melt, fixed density, modes, RPA and the sweep
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

## 6. Running it

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

### eps_AB sweep

```bash
python3 -m melt.epsab_scan --dry-run                                  # 148-run plan
python3 -m melt.epsab_scan --manifest-only --out $OUT --campaign-id epsab_kappa05
python3 -m melt.epsab_scan --run --out $OUT --campaign-id epsab_kappa05 --platform CUDA
python3 -m melt.epsab_scan --status --out $OUT --campaign-id epsab_kappa05
python3 -m melt.epsab_scan --analyze --out $OUT --campaign-id epsab_kappa05
```

`--shard-index K --shard-count S` runs every plan index congruent to `K`
modulo `S` (one shard per scheduler array task); `--parallel P
--threads-per-run T` runs `P` runs at once as subprocesses on a many-core
machine. Both share one campaign directory safely through per-run locks.
`$OUT` must be outside the git tree, since a campaign requires a clean tree.
Stage 2 uses a new `--campaign-id` with a finer `--eps-ABs` list, more
`--seeds` and `--sizes 144 288`. Any run's two-time correlations can also be
inspected directly:

```bash
python3 -m melt.dynamic_structure $OUT/epsab_kappa05/runs/M0144_N40_kappa0p5_pi0p990_epsAB0p5_seed1
python3 scripts/modes_from_trajectory.py --analyze path/to/archived_run   # from a trajectory.npz
```

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

### Sherlock

The `eps_AB` sweep runs on Stanford's Sherlock cluster as a GPU job array.
[`sherlock/README.md`](sherlock/README.md) walks through the SSH setup, the
one-time environment (`python` module plus a pip-installed `openmm[cuda12]`
virtualenv; Sherlock advises against conda), the smoke test, the array
submission, monitoring, and pulling results back through the data transfer
nodes.

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
