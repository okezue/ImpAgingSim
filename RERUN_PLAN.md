# Corrected-Campaign Rerun Plan

Date: 2026-05-22

## What is fixed in the code

| Issue (reviewer item) | File | Status |
|---|---|---|
| Temperature units (#1) | `melt/integrator.py`, `melt/model.py`, `melt/run.py` | T* now converted to Kelvin via `tstar_to_kelvin(T_star, eps_kjmol)`; `LangevinMiddleIntegrator` receives explicit `kelvin` units; unit test asserts T*=0.7 maps to 84.187 K |
| εAB scales both repulsion and attraction (#2) | `melt/integrator.py` | Force field now has shared WCA repulsive core + separately tunable attractive tail. εAB controls only the cross-attraction. Unit tests verify shared core across pair types and εAB-only tail scaling |
| Sequence λ convention mismatch (#3) | `melt/sequences.py` | Single unified `generate_correlated(N, kappa, pi, f_A, rng)` using `P(A→B)=2(1-π)(1-f_A)`, giving λ=2π-1 for every f_A. Empirical λ matches within 0.004 across f_A ∈ {0.3, 0.5, 0.7} |
| Markov sequences across chain boundaries (#3b) | `melt/sequences.py`, `melt/run.py` | New `generate_per_chain` helper generates one independent Markov realisation per chain. Unit test asserts cross-chain-boundary correlation < 0.1 at π=0.95 |
| RPA lag-0 missing (#4) | `paper/sections/01_sequence_control.tex`, `paper/sections/08_methods.tex` | Eq. eq:sequence_autocorr and Eq. eq:rpa_zero corrected to S_0(0) = 4 f_A(1-f_A) [1 + 2κ²λ/(1-λ)]. Methods notes the correction |
| Aging overstatement (#5) | `paper/sections/00_abstract.tex`, `paper/sections/02_aging_inversion.tex` | Abstract no longer claims aging removal. Section retitled "Structural-contrast drift over long waiting times" with absolute-vs-fractional drift framing. Explicit statement that glass-aging claims require F_s, τ_α, MSD, χ_4 which are pending |
| εAB no-ODT overclaim (#5b) | `paper/sections/06_interaction_density.tex` | Removed "no order-disorder boundary" claim; reframed as cross-attraction probe with shared WCA core. ODT claim deferred to finite-size scaling work |

New module `melt/dynamics.py` adds the two-time observables required for a real glass-aging claim: `F_s(k, t; t_w)`, `MSD(t; t_w)`, `α_2(t; t_w)`, self-overlap `Q(t; t_w)`, `χ_4(t)` from realisation variance, and KWW fits.

## How to launch the corrected campaign

1. Push the fix branch to GitHub.
2. From a workstation with AWS credentials and an SSH key:
   ```bash
   export S3_BUCKET=your-bucket
   export GIT_REF=fix-rerun
   SCAN_KIND=rerun ./aws/launch.sh
   ```
   The instance bootstraps the conda env, checks out `${GIT_REF}`, runs `aws/rerun_corrected.sh`, and syncs results to S3.
3. `aws/rerun_corrected.sh` re-runs the decisive panels of Figs 1-5 at corrected reduced units and the corrected force field. The aging campaign now saves per-bead trajectories so `melt.dynamics` can compute F_s, τ_α, MSD, α_2, χ_4 post hoc.

## After the rerun

1. Pull S3 results to `output/aws_corrected/`.
2. Regenerate per-figure CSVs in `paper/source_data/` from the new run directories using the existing analysis scripts (`paper/make_melt_manuscript_assets.py`).
3. Run `python -c "from melt.dynamics import compute_all_dynamics; ..."` per aging chunk to produce F_s(k_*, t; t_w), τ_α(t_w), MSD, χ_4 plots. Add as new Fig 5 panels or as extended data.
4. Update the numerical statements in Sections 2-7 to match the new data. Most prose can stay as the qualitative story is unchanged; numerical values (3.52×, 5.15×, R²=0.83, slopes, etc.) will need to be re-extracted.
5. Re-fit the RPA collapse using the corrected predictor `1 + 2κ² Σ (1-ℓ/N) λ^ℓ` and report the new R² and slope in Fig 2B.
6. Final test pass: `python3 -m pytest tests/test_melt.py -v` should remain green.

## Estimated cost

Per the reviewer's analysis the original campaign was ~83 chunks. The corrected campaign re-runs roughly the same number of conditions but with 5× longer aging trajectories saved for dynamics post-processing. At ~30 minutes/run on an L40S, that is ~250 GPU-hours, ~$1500 at AWS on-demand. Spot instances cut that to ~$500.

## Honest scope

After the rerun the central defensible claim becomes "sequence statistics provide a quantitative design axis for microphase contrast in heteropolymer melts, with κ²/[2(1-π)] as the effective excess-blockiness scale". The aging-inversion headline lives or dies based on F_s, χ_4, and τ_α at multiple waiting times; we will not reinstate that claim without those observables.
