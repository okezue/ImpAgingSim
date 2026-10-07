# README figures

SVGs provide vector artwork and PNGs provide raster copies. Captions are in the
[root README](../../README.md).

## Model and fixed-density artwork

| Asset | Source |
|---|---|
| [Polymer model](polymer-model.svg) | [`render_schematic_figures.py`](../../scripts/render_schematic_figures.py); a schematic 40-bead A/B chain and periodic multichain melt |
| [Sequence construction](sequence-construction.svg) | [`render_schematic_figures.py`](../../scripts/render_schematic_figures.py), [`melt/sequences.py`](../../melt/sequences.py); keep/redraw inheritance and analytic covariance |
| [Interaction potential](interaction-potential.svg) | [`render_potential_figure.py`](../../scripts/render_potential_figure.py), [`melt/integrator.py`](../../melt/integrator.py); the implemented step-gated nonbonded energy |
| [Measured structure](measured-structure.svg) | [`render_structure_figure.py`](../../scripts/render_structure_figure.py); [fixed-density analysis tables](../../output/melt/fixed_density_size/fixed_density_pi099_v2/analysis/) |

The polymer and sequence geometry is schematic. In those drawings, teal denotes
A and orange B. The
mask example includes a redraw that happens to retain the same species. The
analytic covariance law describes the unconditioned generator; conditioning on
an exact global A count can change it.

The potential uses `sigma = eps_core = eps_AA = eps_BB = 1`, `eps_AB = 0.1` and
`r_cut = 2.5 sigma`. Its implemented attractive branch is zero below
`2^(1/6) sigma`, so its onset has an energy jump. Adjacent bonded beads are
excluded from the nonbonded interactions.

The measured-structure figure reads `shell_spectra_condition.csv` and
`condition_summary.csv` from `fixed_density_pi099_v2/analysis`. It displays
`S_psi,psi^(N)/2` with total-bead normalization. Each seed contributes the mean of
its final five saved configurations (steps 242000–250000). Bands and bars are SEM
across five independent seeds. Lines connect measured shell values without
smoothing. The selected peak maximizes the across-seed mean spectrum; it is not
the mean of per-seed maxima. At kappa = 1, the first two system sizes select the
first reciprocal shell and the largest selects the second. The figure does not
extrapolate below the smallest admitted shell or portray a measured morphology.

## Regeneration

From the repository root:

```bash
python scripts/render_schematic_figures.py --export
python scripts/render_potential_figure.py
python scripts/render_structure_figure.py
```

The renderers require NumPy, SciPy and Matplotlib. The schematic `--export`
option also requires Inkscape to outline SVG text and create PNGs. Without
`--export`, the schematic renderer retains editable text. Its optional
`--figma-json-dir DIRECTORY` emits geometry and label coordinates for
reconstructing native Figma text layers.
