# Derived tables and figures

Each study directory below holds the aggregated tables (CSV/JSON), the campaign
`manifest.json` that pins the design and code commit, the figures, and a `README.md`
with the findings. Raw simulation output is not versioned; it is on Zenodo (DOIs
listed at the end of the top-level README).

## The `eps_AB` sweep program (September 2026): transition and dynamics of the melt

Everything at the production geometry (144 x 40 beads, `L = 22 sigma`) with
`eps_AA = eps_BB = 1` and the A-B attraction `eps_AB` swept; `delta_eps = 1 - eps_AB`
is the incompatibility. Zenodo record 10.5281/zenodo.22803471 (version 4.0.0 of the
main archive).

| directory | campaign(s) | runs | question | one-line result |
|---|---|---|---|---|
| [`epsab_stage1/`](epsab_stage1/README.md) | `epsab_kappa05` | 148 | coarse scan of the mixed -> demixed transition at `kappa = 0.5`, `T* = 0.7` | demixing sets in for `0.1 <= delta_eps <= 0.3`; the chi = 0 melt is mixed; the 250k-step protocol does not equilibrate the demixed side |
| [`epsab_stage2/`](epsab_stage2/README.md) | `epsab_stage2_kappa05` | 400 | refine the location, test finite size (144 vs 288 chains) | `eps_AB,c = 0.85 +/- 0.03`, identical in both boxes; RPA spinodal 0.905, `alpha ~ 2`; domain spacing ~14 sigma is not box-limited |
| [`dsf_series/`](dsf_series/README.md) | `dsf_T0p7 ... dsf_T2p0` | 176 | at which temperature does the composition field relax? | `T* >= 1.5` is ergodic (Rouse `q^-4` regime, critical slowing at low q, chi-axis collapse); `T* = 0.7` is arrested |
| [`dsf_theory/`](dsf_theory/README.md) | `dsf_theory_T1p5`, `dsf_theory_T2p0`, `dsf_chi0_M288_T1p5` | 116 | the dataset for comparison with the dynamic theory | single-exponential decay, `tau ~ q^-3.9`, `Gamma(q) S(q) ~ const`; puzzles: relaxation plateau below `q ~ 0.3`, 2.2x speed-up from `T* = 1.5` to 2.0 |
| [`kappa_boundary/`](kappa_boundary/README.md) | `kappa0p0_T1p5 ... kappa1p0_T1p5` | 220 | the phase boundary versus sequence correlation | no transition for `kappa <= 0.25`; `delta_eps_c = 0.20, 0.15, 0.115` for `kappa = 0.5, 0.75, 1`; lamellar signal only at `kappa = 1` |
| [`coarsening/`](coarsening/README.md) | `coarsening_T0p7` | 12 | how the demixed pattern grows over 100,000 tau | `S ~ t^0.5` for three decades near the crossover; deeper quenches arrest at a finite domain size |

## Applications with asymmetric energetics (B loves B, A-B neutral, free volume)

Separate Zenodo records.

| directory | campaign | runs | question | one-line result | DOI |
|---|---|---|---|---|---|
| [`protein_condensation/`](protein_condensation/README.md) | `protein_condensation` | 360 | when does a copolymer solution condense, and into what? | condensation is programmed by sequence blockiness (`kappa >= 0.5`), nearly independent of sticker strength; scattered patterns collapse chain by chain | 10.5281/zenodo.22819768 |
| [`chromatin_memory/`](chromatin_memory/README.md) | `chromatin_memory` | 160 | do transient marked blobs get stabilized by reader-writer feedback? | blobs are transient with a turnover-set lifetime; stabilized once `k_fb / k_off >~ 10`; a stabilized blob then spreads to a fraction set by that ratio | 10.5281/zenodo.22819770 |

Cross-campaign figures and scripts: `scripts/compare_dsf_campaigns.py`,
`scripts/compare_kappa_campaigns.py`, `scripts/coarsening_analysis.py`.

## Paper-era analyses (sections 1-3 of the top-level README)

| directory | contents |
|---|---|
| [`submission_rpa/`](submission_rpa/README.md) | matrix-RPA audit for the Scientific Reports correction, reproduced from the published source tables |
| `figures/` | Fig. S9 finite-size figure (`make_finite_size_figure.py`) |
| `quantitative/`, `comparison/`, `extended/`, `extended_full/` | four-point susceptibility `chi_4`, overlap collapse and ensemble separation tables and figures for the aging campaigns |
| `robustness/`, `robustness_full/` | peak-fit robustness and the kappa scan of `chi_4*` |
| `plots/`, `chi4_*.png` | per-condition `chi_4` curves |
| `quench_evolution/` | scripts for the pre-quench evolution animation |
| `regen/` | regeneration package for the paper figures |
| `summary_condition.csv`, `summary_per_tw.csv` | aging summary tables |

The AWS campaign tables of the paper live in `analysis_aws/`.
