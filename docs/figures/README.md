# README figures

The four figures use a common palette and a white background for legibility in
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
```

The schematic script can also run without `--export` to produce SVGs with live
text. Its optional `--figma-json-dir DIRECTORY` emits geometry and label
coordinates for reconstructing native editable text layers in Figma. No raw
trajectory download or OpenMM installation is needed for these figure scripts.
