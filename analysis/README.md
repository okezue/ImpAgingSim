# Paper analysis and source data

ImpAgingSim contains the analysis supporting **Independent control of
sequence-correlation amplitude tunes post-quench composition fluctuations
in A/B heteropolymer melts**. The corrected manuscript and Supplementary
Information define the reported estimator, normalization, run counts and
interpretation. Follow-up campaigns are maintained in
[HetPoly](https://github.com/okezue/HetPoly).

| Location | Contents and use |
|---|---|
| [`paper/`](paper/README.md) | Final manuscript source and sensitivity tables from the checksum-verified supplementary archive, with provenance and portable figure-regeneration instructions |
| [`submission_rpa/`](submission_rpa/README.md) | Reproduced matrix-RPA audit, formal long-wavelength bare-kernel diagnostic and finite-wavevector predictor fits |
| [`../output/melt/fixed_density_size/fixed_density_pi099_v2/analysis/`](../output/melt/fixed_density_size/fixed_density_pi099_v2/analysis/) | Exact-shell spectra by seed and condition, condition-selected peaks, common-wavevector comparisons and convergence diagnostics for all 30 fixed-density runs |
| [`../docs/corrections/scientific_reports_rpa/`](../docs/corrections/scientific_reports_rpa/) | Corrected manuscript/supplementary sources, figure PDFs, clean and review PDFs, numerical audit and source-table updates |
| [`../docs/figures/`](../docs/figures/README.md) | README schematics, exact exports of manuscript figures and their regeneration instructions |

## Measurement conventions

The balanced primary scans report the peak of
$C_A(k)=S_{AA}^{(N)}(k)-S_{AB}^{(N)}(k)$, with all partial spectra normalized
by total bead count. A run summary averages its final five stored spectra;
independent seeds supply the statistical replicates. At exact A/B exchange
symmetry, $C_A=S_{\psi\psi}^{(N)}/2$; finite samples need not have exactly
equal A and B partial spectra.

The fixed-density comparison measures $S_{\psi\psi}^{(N)}/2$ directly at
nonzero periodic reciprocal vectors and averages within exact shells. Each
condition chooses the maximum of its across-seed mean spectrum, then reports
the seed mean and SEM at that common shell. This differs from averaging
independently selected per-seed maxima. The committed `condition_summary.csv`
records both quantities with explicit column names.

The finite-$k$ predictor fit uses each seed's measured peak-bin centre and
fixed $N=40$, $b=\sigma$. Its $R^2$ describes a response-selected association.
The bare RPA kernel gives a negative formal long-wavelength composition
stiffness; it is a diagnostic of the uncalibrated homogeneous closure and
does not locate the physical melt spinodal. Off-stoichiometric scans require
separation of total-density and composition modes, as detailed in the
Supplementary Information.

## Reproduce the audit and figures

```bash
python3 scripts/paper/figures.py --out output/paper_figures
python3 -m scripts.submission_rpa_audit \
  --source-data analysis/paper/source_data --output /tmp/impaging-rpa
python3 scripts/render_structure_figure.py
```

The [committed final tables](paper/source_data/) reproduce the predictor
fits as well as the bare-kernel diagnostic. Their original source is
`Supplementary_Data_1_source_tables.zip` in the
[paper archive](https://doi.org/10.5281/zenodo.20499120); the ZIP also contains
selected validation trajectories and the original article-wide plotting scripts.

The RPA-only figure regeneration command, source dependencies and manuscript
build commands are documented in the
[correction instructions](../docs/corrections/scientific_reports_rpa/RPA_correction_instructions.md).
The [top-level README](../README.md) gives the production protocol and
checksum-verified fixed-density launcher. Regenerating the static summaries
or the RPA correction requires no new simulation campaign.
