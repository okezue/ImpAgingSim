# README figures

The figures use a common palette and a white background for legibility in
both GitHub themes. SVGs are the canonical README assets; PNGs are convenient
raster exports. Text is outlined in the final SVGs for consistent mathematical
symbols and spacing across browsers. Change labels in the Python sources;
the SVG geometry remains editable in Figma, Inkscape, or Illustrator.

| Figure | What is drawn | Source |
|---|---|---|
| [Polymer model](polymer-model.svg) | A 40-bead chain and periodic continuation in a schematic melt | [`render_schematic_figures.py`](../../scripts/render_schematic_figures.py); geometry follows [`melt/box.py`](../../melt/box.py) and the README's production parameters |
| [Sequence construction](sequence-construction.svg) | Keep/redraw inheritance and analytic normalized covariance | [`render_schematic_figures.py`](../../scripts/render_schematic_figures.py); [`melt/sequences.py`](../../melt/sequences.py) |
| [Interaction potential](interaction-potential.svg) | The exact step-gated nonbonded energy | [`render_potential_figure.py`](../../scripts/render_potential_figure.py); [`melt/integrator.py`](../../melt/integrator.py) |
| [Measured structure](measured-structure.svg) | Exact-shell scattering and peak amplitudes with five-seed SEM | [`render_structure_figure.py`](../../scripts/render_structure_figure.py); [corrected fixed-density campaign tables](../../output/melt/fixed_density_size/fixed_density_pi099_v2/analysis/) |
| [Protein pattern concept](protein-pattern-concept.svg) | Equal-composition chains and an explicit multichain contact topology | [`render_protein_concept.py`](../../scripts/render_protein_concept.py); schematic 40-bead chains, each 12 B / 28 A |
| [Protein condensation results](protein-condensation-results.svg) | Largest-cluster fractions, chain dimensions and local B contacts | [`render_protein_figure.py`](../../scripts/render_protein_figure.py); [protein campaign tables](../../analysis/protein_condensation/) |
| [Chromatin feedback concept](chromatin-feedback-concept.svg) | Mark conversion, spatially local feedback and changing site identity | [`render_chromatin_concept.py`](../../scripts/render_chromatin_concept.py); effective rules in [`melt/marks.py`](../../melt/marks.py) |
| [Chromatin memory results](chromatin-memory-results.svg) | Cluster connectivity, correlation times and mark abundance | [`render_chromatin_figure.py`](../../scripts/render_chromatin_figure.py); [chromatin campaign tables](../../analysis/chromatin_memory/) |

The geometry is schematic, not a trajectory rendering. Teal denotes A and orange
B; purple highlights periodic continuation or a comparison parameter, according
to the local labels. The mask example includes a redraw that retains the same
species, and the covariance plot shows the separate lag-zero value explicitly.
The covariance law describes the unconditioned generator; conditioning on an
exact global A count can change it.

The interaction plot uses `sigma = eps_core = eps_AA = eps_BB = 1`,
`eps_AB = 0.1`, and `r_cut = 2.5 sigma`. The implemented attractive branch is zero
below `2^(1/6) sigma`; its onset therefore has an energy jump. This plot reflects
the existing implementation. Adjacent bonded beads are excluded from these
nonbonded interactions.

The measured figure reads `shell_spectra_condition.csv` and
`condition_summary.csv` from `fixed_density_pi099_v2/analysis`. Its channel is
`S_psi,psi^(N)/2`, with total-bead normalization. Each seed contributes the mean of
its final five saved configurations (steps 242000–250000). Bands and error bars
are the standard error across five independent seeds, not variation across
reciprocal vectors. Lines connect sampled values without smoothing. The selected
peak is the maximum of the across-seed mean spectrum, not the mean of per-seed
maxima. At kappa = 1 the first two system sizes peak at the first reciprocal
shell; the largest peaks at the second shell. No extrapolation below the first
shell or synthetic morphology is shown.

## Regeneration

From the repository root, with NumPy, SciPy, Matplotlib, and Inkscape installed:

```bash
python scripts/render_schematic_figures.py --export
python scripts/render_potential_figure.py
python scripts/render_structure_figure.py
python scripts/render_protein_concept.py
python scripts/render_protein_figure.py
python scripts/render_chromatin_concept.py
python scripts/render_chromatin_figure.py
```

The schematic script can also run without `--export` to produce SVGs with live
text. Its optional `--figma-json-dir DIRECTORY` emits geometry and label
coordinates for reconstructing native editable text layers in Figma. No raw
trajectory download or OpenMM installation is needed for these figure scripts.

## Application figure sources and interpretation

The protein results use 360 runs, grouped into 90 conditions. Each seed contributes
its mean over steps 502,000–1,000,000, with 250 saved cluster summaries. The heatmaps
show the largest B-contact-cluster fraction at both compositions; the lower panels
show radius of gyration and B–B coordination at B fraction 0.3. Color is a condition
mean, with no smoothing or fitted phase boundary. Attraction-axis cell widths follow
the sampled numerical spacing.

The chromatin results use 160 runs, grouped into 40 conditions. The final 6,250 tau
of each run form the analysis window. The two heatmaps use categorical cells at the
sampled turnover and feedback rates (including zero feedback); their shared color
scale is the largest marked-cluster fraction. The correlation-time panel shows
first 1/e times of the B-density correlation at each run's selected spectral peak,
and the time-centered site-mark correlation. The site-mark curve uses eps_BB=1;
with no feedback its kinetics are independent of the spatial neighborhood. The
mark-fraction panel shows every condition separately, without a ratio-only fit;
its horizontal scale is logarithmic above 0.3 and linear near zero.

Both result generators read `per_condition.csv` and independently verify every
plotted mean and SEM against `per_run.csv`. SEM is across four independent seeds,
not across time samples. Contacts include bonded B neighbors within 1.5 sigma.
Density decorrelation is not the survival of a tracked cluster. Persistent density
also does not mean that individual marks remain unchanged, or establish biological
inheritance. All time units are simulation units.

The concept figures are exact, editable vector drawings. They illustrate the
model and do not portray measured morphologies. Protein sequences have verified
bead counts; associated backbones have no crossings or overlapping beads. The
chromatin neighborhood includes contour-distant marked sites. Nucleosome motifs
are visual analogies: the simulation contains spherical polymer beads, with no
explicit DNA, readers, writers or replication. Positions in its paired time
illustrations are intentionally held fixed to isolate changes of site identity.

BioRender was used to explore the [protein concept](https://app.biorender.com/illustrations/87fa61e4275edb336d2cc92a)
and [chromatin concept](https://app.biorender.com/illustrations/9d226c6d7f51ab80995c6b7b).
Those are working drafts. The canonical final figures are the SVGs above, rebuilt
with explicit geometry and outlined labels for reproducible editing and export.

## Zenodo verification

On 2026-09-18, `per_run.csv`, `per_condition.csv`, `summary.json`, `manifest.json`
and the original study README were extracted from both published archives using
validated HTTP byte ranges and checked against commit
`81532887abaec462e712321ada698e6359e99e25`. Every file was byte-identical.
The [verification manifest](application-source-verification.json) records SHA-256
hashes and archive member names. This verifies the extracted members; the full
multi-gigabyte archives were not downloaded to recompute their overall MD5 sums.
The study READMEs have since been revised to distinguish measurements from
interpretation; the numerical source files are unchanged.

- Protein: [Zenodo 22819768, version 1.0.0](https://zenodo.org/records/22819768),
  `protein_condensation_data.tar`.
- Chromatin: [Zenodo 22819770, version 1.0.0](https://zenodo.org/records/22819770),
  `chromatin_memory_data.tar`.

These archives provide cluster summaries, Fourier modes and mark histories;
coordinate trajectories were not present in the sampled completion manifests.
No new figure is presented as a trajectory rendering.
