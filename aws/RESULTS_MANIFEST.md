# Simulation and study-data manifest

The static study includes 606 balanced-composition runs and 168
off-stoichiometric runs used to diagnose density/composition mode mixing.
Three auxiliary trajectories support the point-particle/cloud-in-cell
measurement comparison. Historical scan IDs are retained for compatibility
with the campaign launchers.

| Historical scan or campaign ID | Conditions / runs | Measurement |
|---|---|---|
| `fig1_baseline` | Eight κ values, five seeds: 40 runs | Baseline correlation-amplitude response |
| `fig1_dense` | Eight κ values, four seeds: 32 runs | Density control |
| `fig1_soft` | Eight κ values, four seeds: 32 runs | Interaction-strength control |
| `fig1_short_chain` | Eight κ values, four seeds: 32 runs | Chain-length control |
| `fig2_pi_kappa` | Eight π × six κ values, five seeds: 240 runs | Persistence–amplitude response |
| `fig3_fA_kappa` | Seven compositions × six κ values, four seeds: 168 runs | Density/composition mode-mixing diagnostic |
| `fig4_epsAB` | Four sequence classes × ten cross-attraction values, five seeds: 200 runs | Cross-attraction response |
| `fixed_density_pi099_v2` | 144/288/576 chains × κ=0/1, five seeds: 30 runs | Fixed-density size comparison |

`aws/rerun_corrected.sh` defines the first seven campaigns, and
`aws/fixed_density.sh` defines the size comparison. Across the four κ controls,
π=0.90; the fixed-density series uses π=0.99 and exact global 50:50 composition.
The soft control changes both attractive-tail depth and temperature relative
to the common repulsive core. The current driver uses the core as its
temperature reference; the soft rerun therefore uses `T_equilibrate=2.0` and
`T_quench=0.28` to preserve the study's physical temperatures, including its
33.676 K quench.

## Data and provenance

All 30 completed fixed-density runs, their campaign manifest, run metadata and
per-run hashes are retained at
[`output/melt/fixed_density_size/fixed_density_pi099_v2/`](../output/melt/fixed_density_size/fixed_density_pi099_v2/).
Its `analysis/analysis_manifest.json` records input/output hashes, analysis
provenance and estimator definitions; `analysis/summary.json` records the
completed-run count and condition summaries.

The [Zenodo archive](https://doi.org/10.5281/zenodo.20499120) contains simulation
archives and numerical source tables. The source-table package is pinned to
[version 4.1.0, record 22887076](https://doi.org/10.5281/zenodo.22887076).
Select the simulation families listed above; the archive also contains
material used by subsequent studies.

Static scan summaries use each run's final five stored spectra before averaging
across independent seeds. The fixed-density estimator selects the maximum of
the across-seed mean direct exact-shell `S_ψψ^(N)(q)/2` spectrum, then reports
the five seed values at that common shell. Mesh-normalized and per-bead values
must be converted before combining data across geometries. Peak wavevectors
are limited by the admitted reciprocal shells of each finite box.

## Reanalysis

Reanalyse the completed fixed-density campaign from its retained inputs:

```bash
python -m melt.fixed_density_size_scan --analyze \
  --out output/melt/fixed_density_size --campaign-id fixed_density_pi099_v2
```

Render its measured-structure figure with:

```bash
python scripts/render_structure_figure.py
```
