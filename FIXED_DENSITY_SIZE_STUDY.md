# Fixed-density 144/288/576-chain study

`melt.fixed_density_size_scan` implements the finite-size study requested during review.  Its
default design contains 30 production runs:

- `M = 144, 288, 576` chains;
- `N = 40` beads per chain;
- `L(M) = 22 (M/144)^(1/3)`, giving the same bead density in every system;
- correlated A/B sequences with `pi = 0.99` and `kappa = 0, 1`;
- five independently seeded replicates per condition (`seed = 1, 2, 3, 4, 5`), matching the
  corrected campaign's replicate convention;
- the existing production force field and schedule: `T* = 5.0 -> 0.7`, 30,000
  equilibration steps, 250,000 production steps, and a 2,000-step snapshot interval.

Only the box length and density-analysis grid scale with `M`.  The default grids are 56, 71,
and 89 cells per axis, respectively, keeping the real-space mesh spacing approximately fixed.
All other simulation settings are identical and recorded in the campaign manifest.

The size comparison uses a canonical, exactly balanced sequence ensemble.  For each attempt,
all `M` chain sequences are drawn independently from the original correlated-chain law.  The
complete set is accepted only when it contains exactly `M*N/2` A beads; otherwise the full set
is redrawn, up to a recorded limit of 100,000 attempts.  This is the original disorder law
conditioned only on global composition.  No chain is copied or constructed as another chain's
complement, although independently drawn chains are not required to be unique.  The accepted
attempt count and realized A fraction are stored per run.  Sequence sampling and initial chain
placement use deterministic child streams from `numpy.random.SeedSequence(seed).spawn(2)`, so
different rejection counts do not change the matched placement stream.

## Inspect and initialize

This command prints the complete 30-run plan as JSON and writes nothing.  It does not import or
require OpenMM:

```bash
python -m melt.fixed_density_size_scan --dry-run
```

Create the immutable campaign manifest without starting a simulation:

```bash
python -m melt.fixed_density_size_scan --manifest-only --platform CUDA
```

The default campaign directory is
`output/melt/fixed_density_size/fixed_density_pi099/`.  Reusing a campaign ID with a different
design or simulation configuration is rejected.

Manifest creation and execution require a clean Git tree.  A campaign is pinned to its creation
commit, and cached runs and analysis are rejected if their clean recorded commits differ.

## Run

Run all 30 cases sequentially:

```bash
python -m melt.fixed_density_size_scan --run --platform CUDA
```

For a scheduler array, create the manifest once and assign one zero-based run index to each of
30 tasks:

```bash
python -m melt.fixed_density_size_scan --run --run-index "$SLURM_ARRAY_TASK_ID" --platform CUDA
```

Each run is written to a staging directory.  Only after every required output exists are file
hashes written and the directory atomically promoted to `runs/`.  A repeated command verifies
those hashes and skips completed runs.  If a process was interrupted, first verify that no job
for that run is active, then explicitly archive its partial directory and restart it:

```bash
python -m melt.fixed_density_size_scan --run --run-index 17 \
  --recover-interrupted --platform CUDA
```

The driver restarts interrupted simulations from their beginning; it resumes the campaign at
run boundaries rather than attempting to resume an OpenMM integrator mid-trajectory.

## Structure factors

Every run contains both the existing gridded spectrum (`structure_factor.npz`) and
`direct_structure_factor.npz`.  The latter is computed directly from bead coordinates, without
mesh interpolation, at all periodic-box modes

```text
q = 2*pi*h/L,  h in Z^3,  0 < |q| <= 1.5
```

for ten production frames: five consecutive earlier frames (steps 232,000--240,000 under the
default schedule) followed by the final five frames (steps 242,000--250,000).  The final block
alone defines the estimator; the earlier block is a convergence diagnostic.  The file stores
the individual integer modes and averages over cubic shells of equal `|h|^2`.  The
particle-level normalization is

```text
rho_a(q) = sum_(j in a) exp(i q.r_j)
S_ab(q) = Re[rho_a(q) rho_b(q)*] / N_total
S_psi_psi^(N)(q) = S_CC(q) = S_AA(q) + S_BB(q) - 2 S_AB(q)
primary(q) = S_psi_psi^(N)(q) / 2
```

Because positions and box lengths are in nanometres, all stored `q` and `k*` values are in
`nm^-1`; this unit is also written into the NPZ, CSV, and JSON metadata.

The primary channel is invariant under exchange of the A and B labels, which matters for finite
stochastic systems where `S_AA` and `S_BB` need not be exactly equal in an individual seed.
`S_AA-S_AB` is retained only in explicitly named diagnostic columns.

The dedicated study fixes the direct-mode ceiling at `1.5 nm^-1` so all sizes share the same
physical range.  The CLI rejects another ceiling.

## Aggregate completed runs

After all cases finish, generate deterministic CSV and JSON tables:

```bash
python -m melt.fixed_density_size_scan --analyze
```

The analysis writes:

- `analysis/per_seed_summary.csv`: realized composition, rejection attempts, convergence
  diagnostics, values at the condition-selected shell, and secondary per-seed maxima;
- `analysis/condition_summary.csv`: the primary physical `k*`, amplitude mean and seed SEM at
  that shell, boundary flags, convergence summaries, and secondary sensitivity results;
- `analysis/shell_spectra_by_seed.csv`: early and final per-seed exact-shell spectra and slopes;
- `analysis/shell_spectra_condition.csv`: shell-wise condition means, seed SEMs, and the primary
  selection flag;
- `analysis/common_q_bin_sensitivity.csv`: the a priori sensitivity analysis using the same nine
  physical bins `[0.15,0.30),...,[1.35,1.50] nm^-1` and the same candidate count at every size;
- `analysis/summary.json` and `analysis/analysis_manifest.json`: machine-readable summaries,
  formulas, input hashes, and output hashes.

For each `(M,kappa)` condition, the primary shell is the maximum of the across-seed mean
final-five-frame exact-shell spectrum.  The reported amplitude and uncertainty are then the
mean and SEM across seeds evaluated at that one selected shell.  This avoids defining the main
size result from an average of response-selected per-seed maxima.  Per-seed maxima remain
explicitly secondary.  A selected first or last candidate is flagged as `q_min` or `q_max`.
The common-bin sensitivity repeats the condition-mean selection with identical physical bin
edges and nine candidates for every size.

Convergence fields report the earlier-five to final-five change and the ordinary least-squares
slope across the final five frames, per seed and as condition mean/SEM.  Analysis normally
requires all 30 valid runs.  `--allow-incomplete-analysis` is available for monitoring an active
campaign and records every missing run in `summary.json`.

## Lightweight validation

The focused tests do not execute OpenMM or start production work:

```bash
python -m unittest tests.test_fixed_density_size_scan -v
```
