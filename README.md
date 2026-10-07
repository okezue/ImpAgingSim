# ImpAgingSim

Simulation code and study data for sequence-controlled A/B heteropolymer melts.
The model separates the range of sequence correlations from their amplitude
and measures post-quench composition fluctuations at fixed chemistry, chain
length, density and mean composition. Follow-up studies are maintained in
[HetPoly](https://github.com/okezue/HetPoly).

[Committed run data](output/melt/fixed_density_size/fixed_density_pi099_v2/) ·
[Campaigns and datasets](aws/RESULTS_MANIFEST.md) ·
[Figure sources](docs/figures/README.md) ·
[Data archive](https://doi.org/10.5281/zenodo.20499120)

## Model and sequence construction

The model is a periodic melt of A/B bead–spring chains, simulated with
Langevin dynamics in OpenMM. Harmonic bonds connect adjacent monomers. Every
pair shares a WCA repulsive core; a separate attractive branch gives like
pairs a deeper attraction than A–B pairs. Changing the A–B attraction therefore
preserves excluded volume.

![A bead–spring chain and a periodic multichain melt.](docs/figures/polymer-model.svg)

**Bead–spring geometry.** Teal and orange identify A and B beads. The periodic
boundary translates a chain segment by one box length. This is a schematic;
the baseline production system contains 144 chains of 40 beads.

| Baseline parameter | Value |
|---|---|
| Chains × beads per chain | 144 × 40 = 5,760 beads |
| Box length; bead density | $22\sigma$; $0.54095\sigma^{-3}$ |
| Mean composition | $f_A=0.5$ |
| Attraction depths | $\epsilon_{AA}=\epsilon_{BB}=1$, $\epsilon_{AB}=0.1$ |
| Core depth; attractive cutoff | $\epsilon_{\mathrm{core}}=1$; $2.5\sigma$ |
| Bond stiffness; rest length | $200$; $\sigma$ |
| Equilibration → quench temperature | $T^*=5\rightarrow0.7$ |
| Langevin friction; time step | $1/\tau$; $0.005\tau$ |
| Equilibration; production | 30,000; 250,000 steps |
| Stored spectra; run summary | Every 2,000 steps; mean of final five spectra |

A Markov backbone $z_i$ sets persistence through $\lambda=2\pi-1$. An
independent mask retains each backbone label with probability $\kappa$ and
otherwise redraws it at the same mean composition:
$s_i=b_i z_i+(1-b_i)u_i$. A redraw can produce the same species.

![A Markov backbone, keep/redraw mask and final sequence, with separate effects on correlation amplitude and range.](docs/figures/sequence-construction.svg)

**Independent sequence controls.** For the unconditioned generator,
$\mathrm{Var}(s_i)=4f_A(1-f_A)$, while at positive contour lag
$\mathrm{Cov}(s_i,s_{i+\ell})=4f_A(1-f_A)\kappa^2\lambda^\ell$.
Thus $\pi$ sets the exponential range and $\kappa$ scales its nonzero-lag
amplitude. The plotted covariance curves are analytical illustrations.

Stochastic chains are drawn independently within their own boundaries; they
are neither copies of one sequence nor required to be unique. Their composition
is fixed in expectation. The fixed-density size comparison conditions full
chain sets on an exact global 50:50 composition. Deterministic alternating
and $A_4B_4$ controls repeat the specified sequence on every chain.

## Committed study data

![Measured exact-shell spectra and composition peaks in the matched fixed-density size series.](docs/figures/measured-structure.svg)

**Fixed-density structure.** Left: 144-chain spectra at $\pi=0.99$ and
$T^*=0.7$. Right: the condition-selected peaks across three sizes. The channel
is the total-bead-normalized $S_{\psi\psi}^{(N)}/2$, measured directly at
periodic reciprocal vectors. Bands and bars show mean ± SEM across five seeds.

| Chains $M$ | Box length $L/\sigma$ | $\kappa=0$ peak | $\kappa=1$ peak | Ratio |
|---|---|---|---|---|
| 144 | 22.000 | 6.23 ± 1.22 | 540.9 ± 13.4 | 86.9× |
| 288 | 27.718 | 6.60 ± 0.75 | 773.3 ± 58.6 | 117.2× |
| 576 | 34.923 | 6.83 ± 1.25 | 575.4 ± 39.6 | 84.3× |

Here $L=22(M/144)^{1/3}\sigma$ keeps bead density fixed. For each condition,
the shell maximizing the across-seed mean final-window spectrum is selected,
then all five seed amplitudes are evaluated at that same shell. The correlated
peak lies at $q_{\min}=2\pi/L$ in the two smaller boxes and
$\sqrt2q_{\min}$ in the largest. Its nonmonotonic amplitude and box-dependent
position leave the bulk wavelength and scaling unresolved.

The committed campaign includes its design manifest, metadata, completion
checksums, saved spectra, snapshots and derived condition summaries. The
matrix-RPA implementation is in [`melt/rpa_matrix.py`](melt/rpa_matrix.py);
it returns structure factors only for stable homogeneous kernels. The
baseline bare closure is unstable and does not determine a physical spinodal.

## Run and reproduce

Install Python dependencies and the test runner:

```bash
python3 -m pip install -r requirements.txt pytest
```

Run one baseline production condition (CUDA-capable OpenMM installation
required for `--platform CUDA`; use `CPU` or `Reference` when appropriate):

```bash
python3 -m melt.run \
  --out output/melt/demo --run_id k1 \
  --sequence correlated --kappa 1 --pi 0.99 --f_A 0.5 \
  --n_chains 144 --chain_length 40 --box_size 22 \
  --T_equilibrate 5 --T_quench 0.7 --lj_eps_AB 0.1 \
  --equilibration 30000 --n_steps 250000 --snapshot_interval 2000 \
  --compute_density --grid_size 56 --save_trajectory \
  --seed 1 --platform CUDA
```

`melt.scan` runs the original sequence, persistence–amplitude, composition
and cross-attraction grids; use `python3 -m melt.scan --help` for the grid
arguments and the [campaign manifest](aws/RESULTS_MANIFEST.md) for the production designs.

Plan and launch a fresh 30-run fixed-density campaign from a clean committed
checkout. Store its output outside the checkout so committed study outputs and
their manifests remain intact:

```bash
python3 -m melt.fixed_density_size_scan --dry-run
python3 -m melt.fixed_density_size_scan \
  --manifest-only --out ../study-runs --campaign-id reproduction --platform CUDA
python3 -m melt.fixed_density_size_scan \
  --run --out ../study-runs --campaign-id reproduction --platform CUDA
python3 -m melt.fixed_density_size_scan \
  --analyze --out ../study-runs --campaign-id reproduction
```

The launcher pins code provenance, verifies completed-run checksums and skips
verified runs. `--run-index` selects a scheduler-array task. The committed
completed campaign is
[`fixed_density_pi099_v2`](output/melt/fixed_density_size/fixed_density_pi099_v2/).

Run the code checks and regenerate the README's model and structure figures:

```bash
python3 -m pytest tests/ -q
python3 scripts/render_schematic_figures.py --export
python3 scripts/render_potential_figure.py
python3 scripts/render_structure_figure.py
```

[Figure documentation](docs/figures/README.md) describes the data inputs and
rendering dependencies. [AWS instructions](aws/README.md) cover the static
parameter suite and fixed-density campaign.

## Repository layout

| Location | Contents |
|---|---|
| `melt/` | Simulation engine, sequence generation, observables and RPA implementation |
| `output/melt/fixed_density_size/fixed_density_pi099_v2/` | Completed 30-run campaign and measured summaries |
| `scripts/` | README figure renderers and archive download/upload tools |
| `aws/` | Campaign launchers and dataset manifest |
| `tests/` | Simulation and numerical regression checks |
| `docs/figures/` | README model and measured-structure figures |

## Data availability

Large study datasets are archived at
[10.5281/zenodo.20499120](https://doi.org/10.5281/zenodo.20499120).
The [campaign manifest](aws/RESULTS_MANIFEST.md) identifies the relevant runs.

| Archive file | Study data |
|---|---|
| `heteropolymer_microphase_data.tar` | Production metadata, stored spectra, snapshots and available trajectories |
| `fixed_density_campaign.tar` | Completed 30-run size comparison and earlier execution provenance |
| `Supplementary_Data_1_source_tables.zip` | Numerical source and sensitivity tables, selected validation trajectories and historical analysis scripts |

Download and checksum-verify the simulation archives:

```bash
python3 scripts/fetch_zenodo.py --dest output/zenodo --only \
  heteropolymer_microphase_data.tar fixed_density_campaign.tar
```
