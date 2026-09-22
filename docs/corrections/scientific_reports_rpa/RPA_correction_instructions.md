# Scientific Reports: correction confined to the RPA

No molecular-dynamics rerun is required for this correction. The existing simulation results and all four reported predictor-fit R² values are retained. The numerical audit and updated figure have already been rerun; the commands below reproduce them.

The scalar equation was not algebraically wrong. It was presented without the matrix theory and the symmetry assumptions that justify it. More substantially, the manuscript interpreted an uncalibrated bare interaction integral as a physical stability result. The correction retains RPA as an explicitly approximate homogeneous-reference framework, derives its scalar limit, and removes that unsupported inference. This is not a newly calibrated quantitative prediction for the late post-quench states.

## The corrected theory

With per-bead partial spectra, the balanced Gaussian intrachain matrix is

\[
\Omega(k)=\frac14\begin{pmatrix}D(k)+S_0(k)&D(k)-S_0(k)\\D(k)-S_0(k)&D(k)+S_0(k)\end{pmatrix}.
\]

Here \(D\) is the density form factor and \(S_0\) is the existing sequence-weighted composition form factor. The compressible closure is

\[
\Gamma=\Omega^{-1}+B(k)\mathbf J+\rho\beta\widetilde w(k)\mathbf E,
\qquad \mathbf S_{\rm RPA}=\Gamma^{-1}.
\]

\(\mathbf J\) is the all-ones matrix; \(\mathbf E\) contains the AA, AB and BB attraction depths. \(B(k)\) is an explicit effective repulsive-reference density contribution. It has not been fitted, and the divergent WCA potential is not Fourier transformed. The Gaussian covariance is physically meaningful only when the full inverse-response matrix is positive definite at the retained wavevectors.

For A/B exchange symmetry, the two physical modes decouple:

\[
S_{nn}^{-1}=D^{-1}+B+\tfrac12\rho\beta\widetilde w(\epsilon_{\rm like}+\epsilon_{AB}),
\qquad
S_{\psi\psi}^{-1}=S_0^{-1}-\chi_{\rm bare}/2,
\]

where \(\chi_{\rm bare}=-\rho\beta(\epsilon_{\rm like}-\epsilon_{AB})\widetilde w\). The common repulsion remains in the density mode. Its cancellation from the explicit composition mode follows from symmetry, not dilution. Away from symmetry, the off-diagonal density–composition term must be retained through the Schur complement.

For the plotted contrast, \(C_{\rm sym}=S_{\psi\psi}/2\), so at a common fixed wavevector the symmetric closure would give \(C_{\rm sym}^{-1}=2/S_0-\chi_{\rm bare}\). The existing regression instead fits inverses of condition means with a response-selected wavevector and free slope/intercept. It is an empirical comparison, not that RPA equation or a measurement of χ. Its fitted finite-wavevector slope is 0.1988016913, and its R² is 0.8955978479; neither changes.

The retained bare diagnostic gives χ = 9.499256357 and composition inverse stiffness −3.749628178 in the formal long-wavelength limit. This says the bare homogeneous Gaussian closure fails. It does not establish the actual melt's spinodal, and it is not a measured canonical zero mode. Packing, reference-chain conformations and effective interactions have not been calibrated from a homogeneous melt. That distinction is made consistently in the revised text.

Rumyantsev and Gavrilov use the same matrix-response framework, then integrate a Gaussian fluctuation determinant to study demixing of two distinct, compositionally identical chain populations. This submission concerns monomer-composition spectra within a statistically balanced stochastic sequence ensemble that permits individual-chain composition variation. Their species-demixing law is therefore not transplanted into this study.

## Exactly what to replace

The package contains complete revised sources, so no manual transcription is needed. `RPA_changes.patch` records the exact text changes against the attachment, and `RPA_replacement_text.tex` contains the principal replacement passages for copying into another working version.

1. Replace the manuscript source's `main.tex` and `main_lineno.tex` with the two files in `manuscript_source/`. The review source retains its separate figure-callout and collected-legend layout.
2. Replace `figures/Fig4_sequence_response_theory.pdf` in that source with the corrected copy included there. Also replace the separately uploaded Fig. 4 artwork with `figures/Fig4_sequence_response_theory.pdf`; a PNG and JPEG are provided.
3. Replace the supplementary source's `supplementary.tex` and `supplementary_lineno.tex` with the versions in `supplementary_source/`. The addition derives the matrix response; existing supplementary equation numbers are preserved.
4. In the extracted `Supplementary_Data_1` package, replace only the three matching files supplied under `Supplementary_Data_1_updates/`: `analysis_scripts/make_figures.py`, `analysis_scripts/rebuild_data.py`, and `source_data/Fig6_rpa_q0_diagnostic.csv`. The CSV changes names/provenance, not numerical values. Rezip that same package for submission.
5. Use the four revised PDFs in `pdf/` as the clean/review manuscript and supplementary PDFs, or compile the supplied sources as below. Rezip each updated LaTeX source directory if the journal requires source ZIPs.
6. Adapt the provided `RPA_response_to_reviewer.md` for the response letter.

The changes to the main paper are restricted to the RPA-related abstract sentence, the matrix-response explanation, the reference-chain assumption, the paragraph preceding Eq. (11), Eq. (11) and its interpretation, Fig. 4's RPA caption, the corresponding Discussion and Methods passages, and two pertinent references. Existing simulation protocols, other results, acknowledgements and unrelated artwork are retained. Fig. 4a–c retain their data and content; panel d changes its labels and interpretation only.

## Reproduce the numerical audit

The analysis-only code is in [draft PR #1](https://github.com/okezue/ImpAgingSim/pull/1). From the ImpAgingSim checkout with that branch:

```bash
git fetch origin fix/submission-rpa-clarification
git switch fix/submission-rpa-clarification
python -m scripts.submission_rpa_audit \
  --source-data /absolute/path/to/extracted/Supplementary_Data_1/source_data \
  --output analysis/submission_rpa
python -m pytest -q tests/test_rpa_matrix.py tests/test_rpa.py
```

Use the existing Python environment with NumPy, SciPy, pandas, Matplotlib and pytest. OpenMM and a GPU are not needed for these commands. The audit computes the bare kernel, checks composition stiffness at the formal long-wavelength limit and the first three nonzero reciprocal shells, and reproduces the existing fits. It does not assign an invented density stiffness, generate a density prediction, or continue negative stiffnesses as physical structure factors. The source-data argument is optional if only the kernel/mode diagnostic is wanted.

The repository addition consists of a matrix-RPA module, this focused audit command and focused tests. Existing simulation engines and later chromatin/protein applications are not modified.

## Regenerate Fig. 4 only

From this correction package's directory:

```bash
python Supplementary_Data_1_updates/analysis_scripts/make_figures.py \
  --rpa-only \
  --source-data /absolute/path/to/extracted/Supplementary_Data_1/source_data \
  --output figures
```

Then copy the generated PDF into `manuscript_source/figures/` before recompiling. The original article-wide `rebuild_data.py` relies on historical input directories absent from the submission ZIP; do not run the full data rebuild for this correction. The new `--rpa-only` route reads the already supplied source tables and touches no other figures.

## Compile the PDFs

Run these commands inside `manuscript_source/`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error main_lineno.tex
```

Run these inside `supplementary_source/`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary_lineno.tex
```

Both source packages include their existing class/style files and figure assets. No bibliography download or BibTeX pass is required.

## Sources for the correction

- Rumyantsev, A. M. and Gavrilov, A. A. *Scaling Law for Sequence-Induced Demixing of Compositionally Identical Copolymers*. ACS Macro Letters 15, 589–594 (2026). [doi:10.1021/acsmacrolett.6c00084](https://doi.org/10.1021/acsmacrolett.6c00084); [author manuscript](https://arxiv.org/html/2602.05153v1), especially Eqs. (3), (5)–(8).
- Morse, D. C. and Chung, J. K. *On the chain length dependence of local correlations in polymer melts and a perturbation theory of symmetric polymer blends*. Journal of Chemical Physics 130, 224901 (2009). [doi:10.1063/1.3108460](https://doi.org/10.1063/1.3108460); [author manuscript](https://arxiv.org/abs/0810.5389).
