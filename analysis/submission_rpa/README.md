# Matrix-RPA audit for the Scientific Reports correction

Output of

```bash
python -m scripts.submission_rpa_audit --source-data <Supplementary_Data_1>/source_data --output analysis/submission_rpa
```

run at commit `9ae4621` against the `Supplementary_Data_1_source_tables.zip` published in
Zenodo record 10.5281/zenodo.22803471 (md5 `951454baf7aaa5ebd1329eb9e8efb69c`).

| file | contents |
|---|---|
| `audit.json` | bare interaction kernel `w_tilde(0) = -13.658 sigma^3`, `chi_bare(0) = 9.4993`, composition stiffness at the formal `q -> 0` limit and the first three reciprocal shells, and the reproduced predictor fits (finite-wavevector slope 0.1988016913, `R^2 = 0.8955978479`) |
| `rpa_modes.csv` | the mode-by-mode stiffness table |

These numbers agree with the `audit/` directory of the correction package
(`docs/corrections/scientific_reports_rpa/`) to 2e-15, the floating-point noise of a
different SciPy build. See `docs/corrections/scientific_reports_rpa/RPA_correction_instructions.md`
for what the correction changes and why.
