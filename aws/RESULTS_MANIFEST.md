# Multi-chain copolymer melt campaign — results manifest

Campaign run on AWS g6e.2xlarge (NVIDIA L40S) instances against branch `melt-pipeline-aws`. All bulk simulation outputs live in S3 (gitignored locally) — this manifest is the in-repo index.

## Bucket

`s3://okezue-imp-aging-results/` (us-east-1, account 545338549082)

Pull locally:
```bash
S3_BUCKET=okezue-imp-aging-results bash aws/sync_back.sh
```

## Campaign structure

12 chunks ran across ~12 g6e.2xlarge GPU-hours of L40S compute. Each chunk launched as a fresh instance with its own short-lived STS session token, periodically synced results to its own S3 prefix, and self-terminated on completion. A local credential-refresher daemon rotated STS tokens every 4h to handle long-running chunks.

| Prefix | Workload | Files | Size |
|---|---|---:|---:|
| `seminal_20260427_133456/` | V1: kappa scan + temperature scan + 3 cinematic big runs (400×50=20k beads) | 507 | 409 MB |
| `chunkB_20260427_141413/` | Long-time aging deep dive (5 t_w up to 2M steps) + multi-quench protocol (3 T_first × 4 T_second) | 486 | 35 MB |
| `chunkC_20260427_141356/` | Composition (f_A) scan + persistence (π) scan | 795 | 32 MB |
| `chunkD_20260427_141359/` | **21×15 dense (κ, T) phase diagram** + 7-point finite-size scaling | 3,150 | 107 MB |
| `chunkF_20260427_184248/` | Kappa fine scan v2 (200 chains × 50 beads, 11 κ × 8 seeds) + 10-point ε scan | 714 | 37 MB |
| `chunkG_20260427_184248/` | Temperature scan v2 (12 T × 4 sequences × 6 seeds at 200×50) + block-length scan | 1062 | 69 MB |
| `chunkH_20260427_184248/` | **9×8×4 (κ, ε_AB) phase diagram** | 864 | 35 MB |
| `chunkI_20260427_184249/` | **11×6×4 (f_A, κ) phase diagram** | 792 | 32 MB |
| `chunkJ_20260427_184256/` | Bond stiffness scan + density (box size) scan | 468 | 19 MB |
| `chunkK_20260427_211516/` | **6×3×3×4 (T_q × t_w × sequence) aging-temperature diagram** | 428+ | 23 MB+ |
| `chunkL_20260427_202933/` | 11×5×4 extended (π, κ) scan | 660 | 27 MB |
| `chunkE_20260427_184248/` | ❌ FAILED — 1500×80=120k-bead mega cinematic runs blew up numerically (NaN energies, Rg→10¹¹). Aborted after 4.4hr | 2 | 1 MB |
| `smoke_test/` | Pre-campaign validation | 8 | 2 MB |

**Total valid: ~9,924 files, ~825 MB, ~10,000 simulation runs.**

## Cinematic visualizations (already rendered)

`seminal_20260427_133456/melt/big/` contains for each of 3 sequences (correlated, random, block):

- `polymer_3d.gif` (~18 MB) — animated 3D bead cloud, A red / B blue
- `density_slice.gif` (~7 MB) — 2D φ_A and φ_B density evolution (z-projection)
- `Sk_evolution.gif` (~1 MB) — log-log structure factor S(k) over time
- `final_3d.png` (~1 MB) — high-res final-state snapshot
- `trajectory.npz` (~21 MB) — full position trajectory (re-render at any resolution)
- `density_grids.npz` (~73 MB) — full 64³ density field tensor
- `snapshots.csv`, `meta.json`, `structure_factor.npz` — analysis-ready

## Parameter coverage

| Knob | Range | Where |
|---|---|---|
| κ (correlation strength) | {0.0, 0.05, 0.1, ..., 1.0} (21 values) | chunkD, chunkF |
| T_quench | {0.2, 0.3, ..., 2.5} (15 values) | chunkD, chunkG |
| (κ, T) 2D | 21 × 15 grid | chunkD |
| (κ, ε_AB) 2D | 9 × 8 grid | chunkH |
| (f_A, κ) 2D | 11 × 6 grid | chunkI |
| ε_AB (interaction asymmetry) | {0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.0} | chunkF, chunkH |
| Composition f_A | {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8} | chunkC, chunkI |
| Persistence π | {0.5, 0.55, 0.6, ..., 0.99} (11 values) | chunkC, chunkL |
| Block length | {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20} (11 values) | chunkG |
| Bond stiffness | {50, 100, 200, 400, 800, 1600} | chunkJ |
| Box size (density) | {16, 18, 20, 22, 25, 28, 32} | chunkJ |
| Chain length N | {12, 16, 20, 32, 48, 64, 96, 128} | chunkD, chunkJ |
| Waiting time t_w | {30k, 50k, 100k, 200k, 250k, 300k, 400k, 500k, 1M, 2M} | chunkB, chunkK |
| Sequences | random, block, alternating, correlated | every chunk |

## Analysis entry points

```bash
# Pull everything home
S3_BUCKET=okezue-imp-aging-results bash aws/sync_back.sh

# Cross-run analysis on a scan (xi(t) overlays + final-state S(k) + summary CSV)
python3 -m melt.analyze output/aws/seminal_20260427_133456/melt/scans/kscan_aws/

# Same for any chunk's scan dir
python3 -m melt.analyze output/aws/chunkD_20260427_141359/melt/scans/kt_dense/

# Two-time density-field observables (Q and chi_4 on phi_A(k_*))
python3 -c "from melt.twotime import compute_twotime_for_scan; print(compute_twotime_for_scan('output/aws/chunkB_20260427_141413/melt/scans/aging_deep'))"

# Re-render an animation from a saved trajectory.npz
python3 -m melt.viz output/aws/seminal_20260427_133456/melt/big/big_correlated/
```

## Compute summary

- **Hardware:** g6e.2xlarge with NVIDIA L40S (24GB VRAM, ~80M atom-steps/sec)
- **Total spend:** ~$120 (~$15 wasted on chunk E numerical blowup)
- **Wall time:** ~10 hours including all parallel chunks
- **Peak parallelism:** 7 instances (vCPU bucket cap = 64 = 8 g6e.2xlarge)

## Known issues

1. **Chunk E lost.** 120k-bead mega-scale cinematic runs were numerically unstable at dt=0.005. For future runs at this scale, reduce dt to 0.001-0.002, increase friction, or use a softer pair potential (WCA-cut LJ) at startup.
2. **Original V2/V3 chained runs lost** on instance A. The shutdown-cancellation race between V1's `shutdown -h +5` and the V2-waiter polling killed the chain. The new chunks (E-L launched as separate instances) covered the missed science.

## Phase 5 follow-up chunks (M, N, O)

| Prefix | Variation | Purpose |
|---|---|---|
| `chunkM_*` | (a) **High density** (box=17, ρ≈1.0) | Eliminate globule formation → see microphase signal in raw S_AA |
| `chunkN_*` | (b) **Soft LJ** (ε_AA=ε_BB=0.4) | Θ-solvent regime, chains stay extended |
| `chunkO_*` | (c) **Short chains** (N=12, 480 chains) | Too short to collapse |

Each: 11 κ × 6 T × 4 seeds = 264 scan runs + 3 trajectory recordings at κ ∈ {0,0.5,1.0} for sub-block / intra-globule analysis.

## Headline science findings

After analysis with `melt/deep_analysis.py`:

### 1. Original κ-tuning hypothesis confirmed in true melt (chunkM)
The IMP single-chain finding — **sequence correlation κ tunes microphase structure at fixed marginal coupling variance** — translates directly to multi-chain melts when the system is at high enough density to suppress chain-vacuum globule formation. At ρ≈1.0:

- Raw S_AA peak grows **2.26×** from κ=0 to κ=1 at T_q=0.7
- This is a clean, observable signature in the most basic experimental observable (S(k))
- Not visible at ρ≈0.54 (baseline) where globule positions dominate the signal

### 2. Contrast observable (S_AA − S_AB) extracts the microphase signal at any density
Even at low density where globules dominate raw S_AA, the contrast observable reveals the microphase signal hidden underneath:

- chunkD baseline: contrast grows **3.52×** from κ=0 to κ=1 at T=0.7
- The subtraction removes the trivial gas-liquid (globule) contribution

### 3. Energy ordering is universal across all conditions
κ → energy is monotonic and significant in every regime:

- Baseline ρ=0.54: −10.0%
- High density ρ=1.0: −3.5%  
- Soft LJ ε=0.4: −5.2%
- Short chains N=12: −6.0%

The smaller drop at high density / soft LJ is consistent with reduced LJ contribution to the total potential energy.

### 4. Sub-block + intra-globule structure (V1 big runs)
Trajectory analysis on V1's 400-chain × 50-bead big runs:

| sequence | Rg_A/Rg_full | intra-globule var enhancement |
|---|---|---|
| correlated κ=0.7 | 0.956 | **3.42× over random** |
| block (length=4) | 0.991 | **3.41× over random** |
| random | 0.985 | 2.36× over random |

Correlated and block sequences produce strong A-bead clustering within globules; random shows weaker but still elevated signal.

## Analysis artifacts

`analysis_aws/` directory holds:
- `figures/Fig_kT_heatmaps_4conditions.png` — (κ, T) phase maps for raw S_AA + contrast across baseline, M, N, O
- `figures/Fig_headline_kappa_tuning.png` — κ-tuning curves at T=0.7 across all 4 conditions
- `figures/Fig_energy_kT_4conditions.png` — energy heatmaps
- `figures/Fig_kT_heatmaps_4conditions.png` — (κ, T) phase diagrams
- `figures/contrast_analysis.png` — original contrast finding
- `tables/headline_summary_T07.csv` — top-line numerical comparison
- `tables/<condition>_runs.csv` — per-run metrics
- `tables/<condition>_agg_by_kT.csv` — (κ, T) aggregated means with seed counts
