# Paper figures

These assets document the sequence-controlled heteropolymer melt study. SVGs
provide vector artwork for the repository; PNGs provide raster copies. Captions
belong in the root README and manuscript rather than inside the artwork.

## Manuscript exports

The five `paper-*` figures are unmodified exports of the submitted figure PDFs
in the [corrected manuscript package](../corrections/scientific_reports_rpa/manuscript_source/).
Their axes, data, panel labels and layouts are preserved. The original plotted
peak notation `Ĉ_N` denotes the per-bead estimator called `C_A` in the corrected
manuscript text. Each export is checked against the correction package's source
PDF checksum; [the export manifest](paper-source-verification.json) records the
source path and SHA-256 hashes for the PDF, SVG and PNG.

| Asset | Paper figure | What it shows |
|---|---|---|
| [Sequence statistics](paper-sequences.svg) | [Figure 1](../corrections/scientific_reports_rpa/manuscript_source/figures/Fig1_sequence_construction.pdf) | Representative sequences, within-chain covariance, chain-boundary checks and the finite-chain zero-wavevector form factor |
| [Correlation-amplitude controls](paper-kappa-controls.svg) | [Figure 2](../corrections/scientific_reports_rpa/manuscript_source/figures/Fig2_kappa_controls.pdf) | Peak amplitude, normalized response, selected bin and chain dimensions for the baseline and three controls |
| [Persistence–amplitude response](paper-sequence-response.svg) | [Figure 3](../corrections/scientific_reports_rpa/manuscript_source/figures/Fig3_pi_kappa_response.pdf) | The sampled persistence–amplitude grid, response sections and first-bin selection frequency |
| [Finite-wavevector predictor and RPA diagnostic](paper-rpa-diagnostic.svg) | [Figure 4](../corrections/scientific_reports_rpa/manuscript_source/figures/Fig4_sequence_response_theory.pdf) | The finite-wavevector descriptive regression, its zero-wavevector comparison, form factors and bare-closure stiffness |
| [Cross-attraction response](paper-cross-attraction.svg) | [Figure 5](../corrections/scientific_reports_rpa/manuscript_source/figures/Fig5_cross_attraction.pdf) | Attraction scan, within-seed comparisons and sequence-class endpoints |

The full captions, normalization, run counts and uncertainty definitions are in
[main.tex](../corrections/scientific_reports_rpa/manuscript_source/main.tex).
The empirical peak estimator uses each seed's mean over its final five stored
spectra; error bars are SEM across independent seeds. The finite-wavevector
regression uses the response-selected peak bin and describes the observed data;
it is not a prospective prediction. The negative bare RPA stiffness diagnoses an
uncalibrated homogeneous closure, not a measured spinodal. The cross-attraction
scan does not resolve a nonzero optimum at low attraction.

These exports preserve the supplied manuscript PDFs. The final numerical source
tables are available in [analysis/paper](../../analysis/paper/), with a portable
plotting script for regeneration from those tables. The RPA correction also
retains its updated diagnostic CSV and analysis scripts.

## Model and fixed-density artwork

| Asset | Source |
|---|---|
| [Polymer model](polymer-model.svg) | [`render_schematic_figures.py`](../../scripts/render_schematic_figures.py); a schematic 40-bead A/B chain and periodic multichain melt |
| [Sequence construction](sequence-construction.svg) | [`render_schematic_figures.py`](../../scripts/render_schematic_figures.py), [`melt/sequences.py`](../../melt/sequences.py); keep/redraw inheritance and analytic covariance |
| [Interaction potential](interaction-potential.svg) | [`render_potential_figure.py`](../../scripts/render_potential_figure.py), [`melt/integrator.py`](../../melt/integrator.py); the implemented step-gated nonbonded energy |
| [Measured structure](measured-structure.svg) | [`render_structure_figure.py`](../../scripts/render_structure_figure.py); [corrected fixed-density analysis tables](../../output/melt/fixed_density_size/fixed_density_pi099_v2/analysis/) |

The polymer and sequence geometry is schematic. In those drawings, teal denotes
A and orange B; the manuscript exports retain their own local colour keys. The
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
python scripts/render_paper_figures.py
python scripts/paper/figures.py --out output/paper_figures
python scripts/render_schematic_figures.py --export
python scripts/render_potential_figure.py
python scripts/render_structure_figure.py
```

The manuscript exporter requires Poppler's `pdftocairo`; PNGs are exported at
300 dpi. It validates source PDF hashes and exports SVG glyph outlines. The
portable numerical plotter writes rebuilt paper figures to `output/paper_figures`
without changing these exact manuscript exports. The schematic and measured renderers require NumPy, SciPy and Matplotlib; the
schematic `--export` option also requires Inkscape to outline SVG text and create
PNGs. Without `--export`, the schematic renderer retains editable text. Its
optional `--figma-json-dir DIRECTORY` emits geometry and label coordinates for
reconstructing native Figma text layers. These figure scripts do not require
OpenMM or raw-trajectory downloads.
