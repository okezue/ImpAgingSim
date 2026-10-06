# Paper source tables

The files in `source_data/` are the final published numerical inputs for the five main figures and Supplementary Figures S1–S8. They were copied unchanged from `Supplementary_Data_1_source_tables.zip` in [Zenodo record 22887076, version 4.1.0](https://doi.org/10.5281/zenodo.22887076). The verified archive MD5 is `048a368f761498b245998b6f32cf83fe`; [provenance.json](provenance.json) records the archive and individual file hashes.

Run these commands from the repository root:

```bash
python3 scripts/paper/figures.py --out output/paper_figures
python3 -m scripts.submission_rpa_audit \
  --source-data analysis/paper/source_data --output /tmp/impaging-rpa
```

The figure command writes source artwork to `output/paper_figures/pdf/` and `png/`, and manuscript filenames to `submission/pdf/` and `submission/png/`. It regenerates the five main figures and all eight supplementary figures without changing the supplied tables or manuscript originals. `--source-data` selects another copy of the same source-data folder. Use `--rpa-only` to generate corrected Figure 4 alone.

Most run-level values average the final five saved spectra before averaging across independent seeds. The three auxiliary measurement-validation trajectories contribute their final 25 frames to Supplementary Figure S4; these are distinct from the primary static simulations. Their processed spectra are included here, while the trajectory files remain in the published archive. Frames and reciprocal modes are not independent replicates.

`source_data/fixed_density/` contains the 30-run v2 finite-size source tables, health diagnostics and original aggregate tables. Figure S8 uses the condition-selected exact shell and the common physical interval `0.15 ≤ qσ < 0.30`; all 30 runs are included. Its seed-level means independently reproduce the archived condition means and standard errors.

To recheck those values directly against the retained raw campaign, run:

```bash
python3 scripts/paper/size.py \
  --campaign output/melt/fixed_density_size/fixed_density_pi099_v2 \
  --out /tmp/imp-paper-size/FigS8_finite_size \
  --tables /tmp/imp-paper-size/tables
```

This mode verifies the completion-file checksums and recomputes the final-five-frame quantities. It writes reanalysis tables to the requested output folder. The default invocation of `size.py` instead plots the archived final tables and verifies their aggregation.

The source-data rebuild script is not included: it requires processed inputs outside this repository. Figure regeneration consumes the published final tables directly. Source metadata retains the archive's original provenance labels and cautions; response-selected regressions and the bare RPA diagnostic remain descriptive quantities rather than phase-boundary predictions.

Check the copied tables with `cd analysis/paper && sha256sum -c SHA256SUMS`. `ARCHIVE_SHA256SUMS.txt` preserves the original full-archive manifest for provenance. That manifest has three older entries for files updated in version 4.1.0: the two original analysis scripts and `Fig6_rpa_q0_diagnostic.csv`. The complete archive matches Zenodo's published MD5, every other archive entry matches its SHA-256, and `SHA256SUMS` records the actual bytes of every table copied here.
