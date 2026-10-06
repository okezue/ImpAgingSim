# Paper simulation and source-data manifest

The scope is defined by the manuscript and Supplementary Information in
[`docs/corrections/scientific_reports_rpa/`](../docs/corrections/scientific_reports_rpa/).
Supplementary Table S3 specifies 606 balanced-composition static runs and 168
off-stoichiometric runs used to diagnose density/composition mode mixing.
Three auxiliary trajectories support the point-particle/cloud-in-cell
measurement comparison. Historical scan IDs below are retained for compatibility
with the source-data scripts; their numbers differ from the final figure numbers.

| Historical scan or campaign ID | Conditions / runs | Manuscript use |
|---|---|---|
| `fig1_baseline` | Eight κ values, five seeds: 40 runs | Main Fig. 2; Supplementary Fig. S5 |
| `fig1_dense` | Eight κ values, four seeds: 32 runs | Main Fig. 2; Supplementary Fig. S5 |
| `fig1_soft` | Eight κ values, four seeds: 32 runs | Main Fig. 2; Supplementary Fig. S5 |
| `fig1_short_chain` | Eight κ values, four seeds: 32 runs | Main Fig. 2; Supplementary Fig. S5 |
| `fig2_pi_kappa` | Eight π × six κ values, five seeds: 240 runs | Main Figs. 3–4; Supplementary Figs. S1–S3 |
| `fig3_fA_kappa` | Seven compositions × six κ values, four seeds: 168 runs | Supplementary Fig. S6, mode-mixing diagnostic |
| `fig4_epsAB` | Four sequence classes × ten cross-attraction values, five seeds: 200 runs | Main Fig. 5; Supplementary Fig. S7 |
| `fixed_density_pi099_v2` | 144/288/576 chains × κ=0/1, five seeds: 30 runs | Main finite-size comparison; Supplementary Fig. S8 and Table S4 |

`aws/rerun_corrected.sh` defines the first seven campaigns, and
`aws/fixed_density.sh` defines the size comparison. Across the four κ controls,
π=0.90; the fixed-density series uses π=0.99 and exact global 50:50 composition.
The soft control changes both attractive-tail depth and temperature relative to
the common repulsive core, as specified in Supplementary Table S2. The current
driver uses the core as its temperature reference; the soft rerun therefore uses
`T_equilibrate=2.0` and `T_quench=0.28` to preserve the paper's physical
temperatures, including its 33.676 K quench.

## Archived source and provenance

- The manuscript source, figure PDFs, correction texts and audit outputs are
  bundled in `docs/corrections/scientific_reports_rpa/`, with file checksums in
  `SHA256SUMS.json`.
- All 30 completed fixed-density runs, their campaign manifest, run metadata and
  per-run hashes are retained at
  `output/melt/fixed_density_size/fixed_density_pi099_v2/`. Its
  `analysis/analysis_manifest.json` records input/output hashes, analysis
  provenance and estimator definitions; `analysis/summary.json` records the
  completed-run count and condition summaries.
- The final Supplementary Data figure and sensitivity tables are retained in
  [`analysis/paper/source_data/`](../analysis/paper/source_data/), including the
  fixed-density source tables. The portable figure generator is
  `scripts/paper/figures.py`.
- The manuscript identifies
  [10.5281/zenodo.20499120](https://doi.org/10.5281/zenodo.20499120) as the complete
  simulation archive. The final source package is pinned to
  [version 4.1.0, record 22887076](https://doi.org/10.5281/zenodo.22887076).
  Select the manuscript's simulation families using this manifest; the archive
  also contains material used by subsequent studies.

Static scan summaries use each run's final five stored spectra before averaging
across independent seeds. The primary fixed-density estimator instead selects
the maximum of the across-seed mean direct exact-shell
`S_ψψ^(N)(q)/2` spectrum, then reports the five seed values at that common shell.
Mesh-normalized and per-bead values must be converted before combining data
across geometries. The final manuscript and SI document these conventions and
the finite-box limits on peak-wavevector interpretation.

## Reanalysis

The completed fixed-density campaign can be reanalysed from its retained inputs:

```bash
python -m melt.fixed_density_size_scan --analyze \
  --out output/melt/fixed_density_size --campaign-id fixed_density_pi099_v2
```

The paper figures and corrected RPA diagnostic can be rebuilt from the retained
source tables without OpenMM or new simulations:

```bash
python scripts/paper/figures.py --out output/paper_figures
python -m scripts.submission_rpa_audit --source-data analysis/paper/source_data \
  --output analysis/submission_rpa
```

The correction workflow is documented in
[`RPA_correction_instructions.md`](../docs/corrections/scientific_reports_rpa/RPA_correction_instructions.md).
Its full historical data rebuild requires archived simulation inputs; use the
portable generator above to render the paper figures directly from their source
tables.
